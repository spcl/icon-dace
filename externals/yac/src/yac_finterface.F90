!>
!! @file yac_finterface.F90
!! @brief Implementation of the YAC Fortran user interface
!!
!! @copyright Copyright  (C)  2013 Moritz Hanke <hanke@dkrz.de>
!!                                 Rene Redler <rene.redler@mpimet.mpg.de>
!!
!! @author Moritz Hanke <hanke@dkrz.de>
!!         Rene Redler <rene.redler@mpimet.mpg.de>
!!
!! @par Warranty
!!
!! This code has been tested up to a certain level. Defects and weaknesses,
!! which may be included in the code, do not establish any warranties by the
!! authors. The authors do not make any warranty, express or implied, or assume any
!! liability or responsibility for the use, acquisition or application of this
!! software.
!!
!! @par Fortran User API
!!
!! For a detailed description of the Fortran and C interface please see
!! \ref interface.
!!

!
! Keywords:
! Maintainer: Moritz Hanke <hanke@dkrz.de>
!             Rene Redler <rene.redler@mpimet.mpg.de>
! URL: https://dkrz-sw.gitlab-pages.dkrz.de/yac/
!
! This file is part of YAC.
!
! Redistribution and use in source and binary forms, with or without
! modification, are  permitted provided that the following conditions are
! met:
!
! Redistributions of source code must retain the above copyright notice,
! this list of conditions and the following disclaimer.
!
! Redistributions in binary form must reproduce the above copyright
! notice, this list of conditions and the following disclaimer in the
! documentation and/or other materials provided with the distribution.
!
! Neither the name of the DKRZ GmbH nor the names of its contributors
! may be used to endorse or promote products derived from this software
! without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
! IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
! TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
! PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
! OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
! PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
! LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!


module mo_yac_iso_c_helpers

  public :: yac_internal_cptr2char, yac_check_strlength

  contains

#define YAC_FASSERT(exp, msg) IF (.NOT. exp) call yac_abort_message(msg, __FILE__, __LINE__)

  !> Convertes a C-pointer to a C-string to a Fortran string
  !! @param[in] cptr C-cointer to C-string
  !! @returns Fortran string
  function yac_internal_cptr2char( cptr ) result (string)

    USE, intrinsic :: iso_c_binding, only: c_ptr, c_char, &
        c_f_pointer,c_size_t

    implicit none

    TYPE(c_ptr), intent(in) :: cptr
    CHARACTER(len=:), allocatable :: string
    CHARACTER(kind=c_char), dimension(:), pointer :: chars
    INTEGER(kind=c_size_t) :: i, strlen

    interface
      function strlen_c(str_ptr) bind ( C, name = "strlen" ) result(len)
        use, intrinsic :: iso_c_binding
        type(c_ptr), value      :: str_ptr
        integer(kind=c_size_t)  :: len
      end function strlen_c
    end interface

    strlen = strlen_c(cptr)
    CALL c_f_pointer(cptr, chars, [ strlen ])
    ALLOCATE(character(len=strlen) :: string)
    DO i=1,strlen
      string(i:i) = chars(i)
    END DO
  end function  yac_internal_cptr2char

  !> Checks that a Fortran string does not exceed YAC_MAX_CHARLEN characters
  !! @param[in] in_string Fortran string
  subroutine yac_check_strlength ( in_string )

    use mo_yac_finterface

    implicit none

    character (len=*), intent (in)  :: in_string

    YAC_FASSERT(LEN_TRIM(in_string) < YAC_MAX_CHARLEN, "ERROR(yac_check_strlength): string " // TRIM(in_string) // "exceeds length of YAC_MAX_CHARLEN")

  end subroutine yac_check_strlength

end module mo_yac_iso_c_helpers

module mo_yac_real_to_dble_utils

  public :: send_field_to_dble, &
            send_field_to_dble_single, &
            send_field_to_dble_ptr, &
            recv_field_to_dble, &
            recv_field_to_dble_ptr, &
            recv_field_from_dble, &
            recv_field_from_dble_ptr

contains

  subroutine send_field_to_dble(field_id,        &
                                nbr_hor_points,  &
                                nbr_pointsets,   &
                                collection_size, &
                                send_field,      &
                                send_field_dble, &
                                send_frac_mask,  &
                                send_frac_mask_dble)

    use mo_yac_finterface
    use iso_c_binding, only: c_ptr, c_f_pointer, c_int, c_associated

    implicit none

    interface

      function yac_get_field_put_mask_c2f_c ( field_id ) &
          bind ( c, name='yac_get_field_put_mask_c2f' )

        use, intrinsic :: iso_c_binding, only : c_int, c_ptr

        integer ( kind=c_int ), value :: field_id
        type(c_ptr)                   :: yac_get_field_put_mask_c2f_c

      end function yac_get_field_put_mask_c2f_c

    end interface

    integer, intent (in)  :: field_id
    integer, intent (in)  :: nbr_hor_points
    integer, intent (in)  :: nbr_pointsets
    integer, intent (in)  :: collection_size
    real, intent (in)     :: send_field(nbr_hor_points, &
                                        nbr_pointsets,  &
                                        collection_size)
    double precision, intent (out) :: send_field_dble(nbr_hor_points, &
                                                      nbr_pointsets,  &
                                                      collection_size)
    real, optional, intent (in) :: send_frac_mask(nbr_hor_points, &
                                                  nbr_pointsets,  &
                                                  collection_size)
    double precision, optional, intent (out) :: send_frac_mask_dble(nbr_hor_points, &
                                                                    nbr_pointsets,  &
                                                                    collection_size)

    integer :: i, j, k
    type(c_ptr) :: put_mask_
    type(c_ptr), pointer :: put_mask(:)
    integer(kind=c_int), pointer :: pointset_put_mask(:)

    if (yac_fget_role_from_field_id(field_id) == YAC_EXCHANGE_TYPE_SOURCE) then

      put_mask_ = yac_get_field_put_mask_c2f_c(field_id)
      if (c_associated(put_mask_)) then
        call c_f_pointer(put_mask_, put_mask, shape=[nbr_pointsets])
        do i = 1, collection_size
          do j = 1, nbr_pointsets
            call c_f_pointer( &
                  put_mask(j), pointset_put_mask, shape=[nbr_hor_points])
            do k = 1, nbr_hor_points
              if (pointset_put_mask(k) /= 0) then
                send_field_dble(k, j, i) = dble(send_field(k, j, i))
              else
                send_field_dble(k, j, i) = 0d0
              end if
            end do
          end do
        end do
        if (present(send_frac_mask)) then
          do i = 1, collection_size
            do j = 1, nbr_pointsets
              call c_f_pointer( &
                    put_mask(j), pointset_put_mask, shape=[nbr_hor_points])
              do k = 1, nbr_hor_points
                if (pointset_put_mask(k) /= 0) then
                  send_frac_mask_dble(k, j, i) = dble(send_frac_mask(k, j, i))
                else
                  send_frac_mask_dble(k, j, i) = 0d0
                end if
              end do
            end do
          end do
        end if
      else
        send_field_dble = dble(send_field)
        if (present(send_frac_mask)) then
          send_frac_mask_dble = dble(send_frac_mask)
        end if
      end if
    else
      send_field_dble = 0d0
      if (present(send_frac_mask)) then
        send_frac_mask_dble = 0d0
      end if
    end if
  end subroutine send_field_to_dble

  subroutine send_field_to_dble_single(field_id,        &
                                      nbr_hor_points,  &
                                      collection_size, &
                                      send_field,      &
                                      send_field_dble, &
                                      send_frac_mask,  &
                                      send_frac_mask_dble)

    use mo_yac_finterface
    use iso_c_binding, only: c_ptr, c_f_pointer, c_int, c_associated

    implicit none

    interface

      function yac_get_field_put_mask_c2f_c ( field_id ) &
          bind ( c, name='yac_get_field_put_mask_c2f' )

        use, intrinsic :: iso_c_binding, only : c_int, c_ptr

        integer ( kind=c_int ), value :: field_id
        type(c_ptr)                   :: yac_get_field_put_mask_c2f_c

      end function yac_get_field_put_mask_c2f_c

    end interface

    integer, intent (in)  :: field_id
    integer, intent (in)  :: nbr_hor_points
    integer, intent (in)  :: collection_size
    real, intent (in)     :: send_field(nbr_hor_points, &
                                        collection_size)
    double precision, intent (out) :: send_field_dble(nbr_hor_points, &
                                                      collection_size)
    real, optional, intent (in)     :: send_frac_mask(nbr_hor_points, &
                                                      collection_size)
    double precision, optional, intent (out) :: send_frac_mask_dble(nbr_hor_points, &
                                                                    collection_size)

    integer :: i, j
    type(c_ptr) :: put_mask_
    type(c_ptr), pointer :: put_mask(:)
    integer(kind=c_int), pointer :: pointset_put_mask(:)

    if (yac_fget_role_from_field_id(field_id) == YAC_EXCHANGE_TYPE_SOURCE) then

      put_mask_ = yac_get_field_put_mask_c2f_c(field_id)
      if (c_associated(put_mask_)) then
        call c_f_pointer(put_mask_, put_mask, shape=[1])
        do i = 1, collection_size
          call c_f_pointer( &
                  put_mask(1), pointset_put_mask, shape=[nbr_hor_points])
          do j = 1, nbr_hor_points
            if (pointset_put_mask(j) /= 0) then
              send_field_dble(j, i) = dble(send_field(j, i))
            else
              send_field_dble(j, i) = 0d0
            end if
          end do
        end do
        if (present(send_frac_mask)) then
          do i = 1, collection_size
            call c_f_pointer( &
                    put_mask(1), pointset_put_mask, shape=[nbr_hor_points])
            do j = 1, nbr_hor_points
              if (pointset_put_mask(j) /= 0) then
                send_frac_mask_dble(j, i) = dble(send_frac_mask(j, i))
              else
                send_frac_mask_dble(j, i) = 0d0
              end if
            end do
          end do
        end if
      else
        send_field_dble = dble(send_field)
        if (present(send_frac_mask)) then
          send_frac_mask_dble = dble(send_frac_mask)
        end if
      end if
    else
      send_field_dble = 0d0
      if (present(send_frac_mask)) then
        send_frac_mask_dble = 0d0
      end if
    end if
  end subroutine send_field_to_dble_single

  subroutine send_field_to_dble_ptr(field_id,        &
                                    nbr_pointsets,   &
                                    collection_size, &
                                    send_field,      &
                                    send_field_dble, &
                                    send_frac_mask,  &
                                    send_frac_mask_dble)

    use mo_yac_finterface
    use iso_c_binding, only: c_ptr, c_f_pointer, c_int, c_associated

    implicit none

    interface

      function yac_get_field_put_mask_c2f_c ( field_id ) &
          bind ( c, name='yac_get_field_put_mask_c2f' )

        use, intrinsic :: iso_c_binding, only : c_int, c_ptr

        integer ( kind=c_int ), value :: field_id
        type(c_ptr)                   :: yac_get_field_put_mask_c2f_c

      end function yac_get_field_put_mask_c2f_c

    end interface

    integer, intent (in)            :: field_id
    integer, intent (in)            :: nbr_pointsets
    integer, intent (in)            :: collection_size
    type(yac_real_ptr), intent (in) :: send_field(nbr_pointsets, &
                                                  collection_size)
    type(yac_dble_ptr), intent (out) :: send_field_dble(nbr_pointsets, &
                                                        collection_size)
    type(yac_real_ptr), optional, intent (in) :: send_frac_mask(nbr_pointsets, &
                                                                collection_size)
    type(yac_dble_ptr), optional, intent (out) :: send_frac_mask_dble(nbr_pointsets, &
                                                                      collection_size)

    integer :: i, j, k, nbr_hor_points
    type(c_ptr) :: put_mask_
    type(c_ptr), pointer :: put_mask(:)
    integer(kind=c_int), pointer :: pointset_put_mask(:)

    if (yac_fget_role_from_field_id(field_id) == YAC_EXCHANGE_TYPE_SOURCE) then

      put_mask_ = yac_get_field_put_mask_c2f_c(field_id)
      if (c_associated(put_mask_)) then
        call c_f_pointer(put_mask_, put_mask, shape=[nbr_pointsets])
        do i = 1, collection_size
          do j = 1, nbr_pointsets
            nbr_hor_points = size(send_field(j,i)%p)
            allocate(send_field_dble(j,i)%p(nbr_hor_points))
            call c_f_pointer( &
                  put_mask(j), pointset_put_mask, shape=[nbr_hor_points])
            do k = 1, nbr_hor_points
              if (pointset_put_mask(k) /= 0) then
                send_field_dble(j, i)%p(k) = dble(send_field(j, i)%p(k))
              else
                send_field_dble(j, i)%p(k) = 0d0
              end if
            end do
          end do
        end do
        if (present(send_frac_mask)) then
          do i = 1, collection_size
            do j = 1, nbr_pointsets
              nbr_hor_points = size(send_frac_mask(j,i)%p)
              allocate(send_frac_mask_dble(j,i)%p(nbr_hor_points))
              call c_f_pointer( &
                    put_mask(j), pointset_put_mask, shape=[nbr_hor_points])
              do k = 1, nbr_hor_points
                if (pointset_put_mask(k) /= 0) then
                  send_frac_mask_dble(j, i)%p(k) = dble(send_frac_mask(j, i)%p(k))
                else
                  send_frac_mask_dble(j, i)%p(k) = 0d0
                end if
              end do
            end do
          end do
        end if
      else
        do i = 1, collection_size
          do j = 1, nbr_pointsets
            nbr_hor_points = size(send_field(j,i)%p)
            allocate(send_field_dble(j,i)%p(nbr_hor_points))
            send_field_dble(j,i)%p = dble(send_field(j,i)%p)
          end do
        end do
        if (present(send_frac_mask)) then
          do i = 1, collection_size
            do j = 1, nbr_pointsets
              nbr_hor_points = size(send_frac_mask(j,i)%p)
              allocate(send_frac_mask_dble(j,i)%p(nbr_hor_points))
              send_frac_mask_dble(j,i)%p = dble(send_frac_mask(j,i)%p)
            end do
          end do
        end if
      end if
    else
      do i = 1, collection_size
        do j = 1, nbr_pointsets
          nbr_hor_points = size(send_field(j,i)%p)
          allocate(send_field_dble(j,i)%p(nbr_hor_points))
          send_field_dble(j,i)%p = 0d0
        end do
      end do
      if (present(send_frac_mask)) then
        do i = 1, collection_size
          do j = 1, nbr_pointsets
            nbr_hor_points = size(send_frac_mask(j,i)%p)
            allocate(send_frac_mask_dble(j,i)%p(nbr_hor_points))
            send_frac_mask_dble(j,i)%p = 0d0
          end do
        end do
      end if
    end if
  end subroutine send_field_to_dble_ptr

  ! -----------------------------------------------------------------------

  subroutine recv_field_to_dble(field_id,        &
                                nbr_hor_points,  &
                                collection_size, &
                                recv_field,      &
                                recv_field_dble)

    use mo_yac_finterface
    use iso_c_binding, only: c_ptr, c_f_pointer, c_int, c_associated

    implicit none

    interface

      function yac_get_field_get_mask_c2f_c ( field_id ) &
          bind ( c, name='yac_get_field_get_mask_c2f' )

        use, intrinsic :: iso_c_binding, only : c_int, c_ptr

        integer ( kind=c_int ), value :: field_id
        type(c_ptr)                   :: yac_get_field_get_mask_c2f_c

      end function yac_get_field_get_mask_c2f_c

    end interface

    integer, intent (in)  :: field_id
    integer, intent (in)  :: nbr_hor_points
    integer, intent (in)  :: collection_size
    real, intent (in)     :: recv_field(nbr_hor_points, &
                                        collection_size)
    double precision, intent (out) :: recv_field_dble(nbr_hor_points, &
                                                      collection_size)

    integer :: i, j
    type(c_ptr) :: get_mask_
    integer(kind=c_int), pointer :: get_mask(:)

    if (yac_fget_role_from_field_id(field_id) == YAC_EXCHANGE_TYPE_TARGET) then

      get_mask_ = yac_get_field_get_mask_c2f_c(field_id)
      if (c_associated(get_mask_)) then
        call c_f_pointer(get_mask_, get_mask, shape=[nbr_hor_points])
        do i = 1, collection_size
          do j = 1, nbr_hor_points
            if (get_mask(j) /= 0) then
              recv_field_dble(j, i) = dble(recv_field(j, i))
            else
              recv_field_dble(j, i) = 0d0
            end if
          end do
        end do
      else
        recv_field_dble = dble(recv_field)
      end if
    else
      recv_field_dble = 0d0
    end if
  end subroutine recv_field_to_dble

  subroutine recv_field_to_dble_ptr(field_id,        &
                                    collection_size, &
                                    recv_field,      &
                                    recv_field_dble)

    use mo_yac_finterface
    use iso_c_binding, only: c_ptr, c_f_pointer, c_int, c_associated

    implicit none

    interface

      function yac_get_field_get_mask_c2f_c ( field_id ) &
          bind ( c, name='yac_get_field_get_mask_c2f' )

        use, intrinsic :: iso_c_binding, only : c_int, c_ptr

        integer ( kind=c_int ), value :: field_id
        type(c_ptr)                   :: yac_get_field_get_mask_c2f_c

      end function yac_get_field_get_mask_c2f_c

    end interface

    integer, intent (in)             :: field_id
    integer, intent (in)             :: collection_size
    type(yac_real_ptr), intent (in)  :: recv_field(collection_size)
    type(yac_dble_ptr), intent (out) :: recv_field_dble(collection_size)

    integer :: i, j, nbr_hor_points
    type(c_ptr) :: get_mask_
    integer(kind=c_int), pointer :: get_mask(:)

    if (yac_fget_role_from_field_id(field_id) == YAC_EXCHANGE_TYPE_TARGET) then

      get_mask_ = yac_get_field_get_mask_c2f_c(field_id)
      if (c_associated(get_mask_) .and. (collection_size > 0)) then
        nbr_hor_points = size(recv_field(1)%p)
        call c_f_pointer(get_mask_, get_mask, shape=[nbr_hor_points])
        do i = 1, collection_size
          nbr_hor_points = size(recv_field(i)%p)
          allocate(recv_field_dble(i)%p(nbr_hor_points))
          do j = 1, nbr_hor_points
            if (get_mask(j) /= 0) then
              recv_field_dble(i)%p(j) = dble(recv_field(i)%p(j))
            else
              recv_field_dble(i)%p(j) = 0d0
            end if
          end do
        end do
      else
        do i = 1, collection_size
          nbr_hor_points = size(recv_field(i)%p)
          allocate(recv_field_dble(i)%p(nbr_hor_points))
          recv_field_dble(i)%p = dble(recv_field(i)%p)
        end do
      end if
    else
      do i = 1, collection_size
        nbr_hor_points = size(recv_field(i)%p)
        allocate(recv_field_dble(i)%p(nbr_hor_points))
        recv_field_dble(i)%p = 0d0
      end do
    end if
  end subroutine recv_field_to_dble_ptr

  ! -----------------------------------------------------------------------

  subroutine recv_field_from_dble(field_id,        &
                                  nbr_hor_points,  &
                                  collection_size, &
                                  recv_field_dble, &
                                  recv_field)

    use mo_yac_finterface
    use iso_c_binding, only: c_ptr, c_f_pointer, c_int, c_associated

    implicit none

    interface

      function yac_get_field_get_mask_c2f_c ( field_id ) &
          bind ( c, name='yac_get_field_get_mask_c2f' )

        use, intrinsic :: iso_c_binding, only : c_int, c_ptr

        integer ( kind=c_int ), value :: field_id
        type(c_ptr)                   :: yac_get_field_get_mask_c2f_c

      end function yac_get_field_get_mask_c2f_c

    end interface

    integer, intent (in)  :: field_id
    integer, intent (in)  :: nbr_hor_points
    integer, intent (in)  :: collection_size
    double precision, intent (in) :: recv_field_dble(nbr_hor_points, &
                                                    collection_size)
    real, intent (inout)    :: recv_field(nbr_hor_points, &
                                          collection_size)

    integer :: i, j
    type(c_ptr) :: get_mask_
    integer(kind=c_int), pointer :: get_mask(:)

    if (yac_fget_role_from_field_id(field_id) == YAC_EXCHANGE_TYPE_TARGET) then

      get_mask_ = yac_get_field_get_mask_c2f_c(field_id)
      if (c_associated(get_mask_)) then
        call c_f_pointer(get_mask_, get_mask, shape=[nbr_hor_points])
        do i = 1, collection_size
          do j = 1, nbr_hor_points
            if (get_mask(j) /= 0) then
              recv_field(j, i) = real(recv_field_dble(j, i))
            end if
          end do
        end do
      else
        recv_field = real(recv_field_dble)
      end if
    end if
  end subroutine recv_field_from_dble

  subroutine recv_field_from_dble_ptr(field_id,        &
                                      collection_size, &
                                      recv_field_dble, &
                                      recv_field)

    use mo_yac_finterface
    use iso_c_binding, only: c_ptr, c_f_pointer, c_int, c_associated

    implicit none

    interface

      function yac_get_field_get_mask_c2f_c ( field_id ) &
          bind ( c, name='yac_get_field_get_mask_c2f' )

        use, intrinsic :: iso_c_binding, only : c_int, c_ptr

        integer ( kind=c_int ), value :: field_id
        type(c_ptr)                   :: yac_get_field_get_mask_c2f_c

      end function yac_get_field_get_mask_c2f_c

    end interface

    integer, intent (in)               :: field_id
    integer, intent (in)               :: collection_size
    type(yac_dble_ptr), intent (inout) :: recv_field_dble(collection_size)
    type(yac_real_ptr), intent (inout) :: recv_field(collection_size)

    integer :: i, j, nbr_hor_points
    type(c_ptr) :: get_mask_
    integer(kind=c_int), pointer :: get_mask(:)

    if (yac_fget_role_from_field_id(field_id) == YAC_EXCHANGE_TYPE_TARGET) then

      get_mask_ = yac_get_field_get_mask_c2f_c(field_id)
      if (c_associated(get_mask_) .and. (collection_size > 0)) then
        nbr_hor_points = size(recv_field(1)%p)
        call c_f_pointer(get_mask_, get_mask, shape=[nbr_hor_points])
        do i = 1, collection_size
          nbr_hor_points = size(recv_field(i)%p)
          do j = 1, nbr_hor_points
            if (get_mask(j) /= 0) then
              recv_field(i)%p(j) = real(recv_field_dble(i)%p(j))
            end if
          end do
          deallocate(recv_field_dble(i)%p)
        end do
      else
        do i = 1, collection_size
          recv_field(i)%p = real(recv_field_dble(i)%p)
          deallocate(recv_field_dble(i)%p)
        end do
      end if
    end if
  end subroutine recv_field_from_dble_ptr

end module mo_yac_real_to_dble_utils

! -------------------------- init -------------------------------------

subroutine yac_fmpi_handshake ( comm, group_names, group_comms )
  use, intrinsic :: iso_c_binding, only : c_ptr, c_null_char, c_loc
  use mo_yac_finterface, dummy => yac_fmpi_handshake

  implicit none

  interface
     subroutine yac_cmpi_handshake_c (comm, n, group_names, group_comms) &
          bind ( c, name='yac_cmpi_handshake_f2c')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       use mo_yac_finterface, only : YAC_MPI_FINT_KIND
       integer (kind = YAC_MPI_FINT_KIND ), intent(in), value :: comm
       integer(c_int), intent(in), value :: n
       type (c_ptr) , intent(in) :: group_names(n)
       integer (kind = YAC_MPI_FINT_KIND ), intent(out) :: group_comms(n)
     end subroutine yac_cmpi_handshake_c
  end interface

  integer, intent(in) :: comm
  character(len=YAC_MAX_CHARLEN), intent(in) :: group_names(:)
  integer, intent(out) :: group_comms(SIZE(group_names))

  CHARACTER (kind=c_char, len=YAC_MAX_CHARLEN+1), TARGET :: &
    group_names_cpy(SIZE(group_names))
  type( c_ptr ) :: group_name_ptr(SIZE(group_names))
  integer :: i
  DO i=1,SIZE(group_names)
     group_names_cpy(i) = TRIM(group_names(i)) // c_null_char
     group_name_ptr(i) = c_loc(group_names_cpy(i))
  END DO

  call yac_cmpi_handshake_c( &
    comm, SIZE(group_names), group_name_ptr, group_comms)

end subroutine yac_fmpi_handshake

subroutine yac_finit_emitter_flags()

  use mo_yac_finterface, only : YAC_YAML_EMITTER_DEFAULT_F, &
                                YAC_YAML_EMITTER_JSON_F

  implicit none

  interface

      function yac_cyaml_get_emitter_flag_default () &
        bind ( c, name='yac_cyaml_get_emitter_flag_default_c2f' )

        use, intrinsic :: iso_c_binding, only : c_int

        integer ( kind=c_int) :: yac_cyaml_get_emitter_flag_default

      end function yac_cyaml_get_emitter_flag_default

      function yac_cyaml_get_emitter_flag_json () &
        bind ( c, name='yac_cyaml_get_emitter_flag_json_c2f' )

        use, intrinsic :: iso_c_binding, only : c_int

        integer ( kind=c_int) :: yac_cyaml_get_emitter_flag_json

      end function yac_cyaml_get_emitter_flag_json

  end interface

  YAC_YAML_EMITTER_DEFAULT_F = yac_cyaml_get_emitter_flag_default()
  YAC_YAML_EMITTER_JSON_F = yac_cyaml_get_emitter_flag_json()

end subroutine

subroutine yac_finit_comm ( comm )

  use, intrinsic :: iso_c_binding, only : c_null_char
  use mo_yac_finterface, dummy => yac_finit_comm

  implicit none

  interface

      subroutine yac_cinit_comm_c ( comm ) &
        bind ( c, name='yac_cinit_comm_f2c' )

        use, intrinsic :: iso_c_binding, only : c_char
        use mo_yac_finterface, only : YAC_MPI_FINT_KIND

        integer ( kind=YAC_MPI_FINT_KIND ), value :: comm

      end subroutine yac_cinit_comm_c

      subroutine yac_finit_emitter_flags()
      end subroutine yac_finit_emitter_flags

  end interface

  integer, intent(in) :: comm

  call yac_finit_emitter_flags ( )
  call yac_cinit_comm_c ( comm )

end subroutine yac_finit_comm

subroutine yac_finit_comm_instance( comm, yac_instance_id)

  use, intrinsic :: iso_c_binding, only : c_null_char
  use mo_yac_finterface, dummy => yac_finit_comm_instance

  implicit none

  interface

      subroutine yac_cinit_comm_instance_c ( comm, yac_instance_id) &
        bind ( c, name='yac_cinit_comm_instance_f2c' )

        use, intrinsic :: iso_c_binding, only : c_char, c_int
        use mo_yac_finterface, only : YAC_MPI_FINT_KIND

        integer ( kind=YAC_MPI_FINT_KIND ), value :: comm
        integer (kind=c_int)                      :: yac_instance_id

      end subroutine yac_cinit_comm_instance_c

      subroutine yac_finit_emitter_flags()
      end subroutine yac_finit_emitter_flags

  end interface

  integer, intent(in)  :: comm
  integer, intent(out) :: yac_instance_id

  call yac_finit_emitter_flags ( )
  call yac_cinit_comm_instance_c ( comm, yac_instance_id )

end subroutine yac_finit_comm_instance

subroutine yac_finit (  )

  use, intrinsic :: iso_c_binding, only : c_null_char
  use mo_yac_finterface, dummy => yac_finit

  implicit none

  interface

    subroutine yac_cinit_c (  ) &
      bind ( c, name='yac_cinit' )

      use, intrinsic :: iso_c_binding, only : c_char

    end subroutine yac_cinit_c

    subroutine yac_finit_emitter_flags()
    end subroutine yac_finit_emitter_flags

  end interface

  call yac_finit_emitter_flags ( )
  call yac_cinit_c ( )

end subroutine yac_finit

subroutine yac_finit_instance ( yac_instance_id )

  use, intrinsic :: iso_c_binding, only : c_null_char
  use mo_yac_finterface, dummy => yac_finit_instance

  implicit none

  interface

      subroutine yac_cinit_instance_c ( yac_instance_id ) &
        bind ( c, name='yac_cinit_instance' )

        use, intrinsic :: iso_c_binding, only : c_char, c_int

        integer (kind=c_int)                  :: yac_instance_id

      end subroutine yac_cinit_instance_c

      subroutine yac_finit_emitter_flags()
      end subroutine yac_finit_emitter_flags

  end interface

  integer, intent(out)         :: yac_instance_id !< [OUT] returned handle to the YAC instance

  call yac_finit_emitter_flags ( )
  call yac_cinit_instance_c ( yac_instance_id )

end subroutine yac_finit_instance

subroutine yac_finit_comm_dummy ( world_comm )

  use mo_yac_finterface, dummy => yac_finit_comm_dummy

  implicit none

  interface

      subroutine yac_cinit_comm_dummy_c ( world_comm ) &
        bind ( c, name='yac_cinit_comm_dummy_f2c' )

        use mo_yac_finterface, only : YAC_MPI_FINT_KIND

        integer ( kind=YAC_MPI_FINT_KIND ), value :: world_comm

      end subroutine yac_cinit_comm_dummy_c

      subroutine yac_finit_emitter_flags()
      end subroutine yac_finit_emitter_flags

  end interface

  integer, intent(in) :: world_comm !< [IN] MPI world communicator (optional)

  call yac_finit_emitter_flags ( )
  call yac_cinit_comm_dummy_c ( world_comm )

end subroutine yac_finit_comm_dummy

subroutine yac_finit_dummy ( )

  use mo_yac_finterface, dummy => yac_finit_dummy

  implicit none

  interface

      subroutine yac_cinit_dummy_c ( ) &
        bind ( c, name='yac_cinit_dummy' )

      end subroutine yac_cinit_dummy_c

      subroutine yac_finit_emitter_flags()
      end subroutine yac_finit_emitter_flags

  end interface

  call yac_finit_emitter_flags ( )
  call yac_cinit_dummy_c ( )

end subroutine yac_finit_dummy


! -------------------------- reading config file ---------------------------

subroutine yac_fread_config_yaml (yaml_filename)
  use mo_yac_finterface, dummy => yac_fread_config_yaml
  use, intrinsic :: iso_c_binding, only : c_null_char
  use mo_yac_iso_c_helpers

  implicit none

  interface
     subroutine yac_cread_config_yaml_c(yaml_filename) &
          bind ( c, name='yac_cread_config_yaml' )
       use, intrinsic :: iso_c_binding, only : c_char
       character (kind=c_char), dimension(*) :: yaml_filename
     end subroutine yac_cread_config_yaml_c
  end interface

  character(len=*), intent(in) :: yaml_filename

  call yac_cread_config_yaml_c(TRIM(yaml_filename) // c_null_char)

end subroutine yac_fread_config_yaml

subroutine yac_fread_config_yaml_instance(yac_instance_id, yaml_filename)
  use mo_yac_finterface, dummy => yac_fread_config_yaml_instance
  use, intrinsic :: iso_c_binding, only : c_null_char
  use mo_yac_iso_c_helpers

  implicit none

  interface
     subroutine yac_cread_config_yaml_instance_c( &
          yac_instance_id, yaml_filename) &
          bind ( c, name='yac_cread_config_yaml_instance' )
       use, intrinsic :: iso_c_binding, only : c_char, c_int
       integer (kind=c_int), value           :: yac_instance_id
       character (kind=c_char), dimension(*) :: yaml_filename
     end subroutine yac_cread_config_yaml_instance_c
  end interface

  integer, intent(in)          :: yac_instance_id
  character(len=*), intent(in) :: yaml_filename

  call yac_cread_config_yaml_instance_c(yac_instance_id, &
       & TRIM(yaml_filename) // c_null_char)

end subroutine yac_fread_config_yaml_instance

subroutine yac_fread_config_json (json_filename)
  use mo_yac_finterface, dummy => yac_fread_config_json
  use, intrinsic :: iso_c_binding, only : c_null_char
  use mo_yac_iso_c_helpers

  implicit none

  interface
     subroutine yac_cread_config_json_c(json_filename) &
          bind ( c, name='yac_cread_config_json' )
       use, intrinsic :: iso_c_binding, only : c_char
       character (kind=c_char), dimension(*) :: json_filename
     end subroutine yac_cread_config_json_c
  end interface

  character(len=*), intent(in) :: json_filename

  call yac_cread_config_json_c(TRIM(json_filename) // c_null_char)

end subroutine yac_fread_config_json

subroutine yac_fread_config_json_instance(yac_instance_id, json_filename)
  use mo_yac_finterface, dummy => yac_fread_config_json_instance
  use, intrinsic :: iso_c_binding, only : c_null_char
  use mo_yac_iso_c_helpers

  implicit none

  interface
     subroutine yac_cread_config_json_instance_c( &
          yac_instance_id, json_filename) &
          bind ( c, name='yac_cread_config_json_instance' )
       use, intrinsic :: iso_c_binding, only : c_char, c_int
       integer (kind=c_int), value           :: yac_instance_id
       character (kind=c_char), dimension(*) :: json_filename
     end subroutine yac_cread_config_json_instance_c
  end interface

  integer, intent(in)          :: yac_instance_id
  character(len=*), intent(in) :: json_filename

  call yac_cread_config_json_instance_c(yac_instance_id, &
       & TRIM(json_filename) // c_null_char)

end subroutine yac_fread_config_json_instance

! -------------------------- cleanup -----------------------------------

subroutine yac_fcleanup ( )

  use mo_yac_finterface, dummy => yac_fcleanup

  implicit none

  interface

      subroutine yac_ccleanup_c () bind ( c, name='yac_ccleanup' )
      end subroutine yac_ccleanup_c

  end interface

  call yac_ccleanup_c ( )

end subroutine yac_fcleanup

subroutine yac_fcleanup_instance ( yac_instance_id )

  use mo_yac_finterface, dummy => yac_fcleanup_instance

  implicit none

  interface

      subroutine yac_ccleanup_instance_c ( yac_instance_id ) &
        bind ( c, name='yac_ccleanup_instance' )

        use, intrinsic :: iso_c_binding, only : c_int

        integer (kind=c_int), value :: yac_instance_id

      end subroutine yac_ccleanup_instance_c

  end interface

  integer, intent(in) :: yac_instance_id !< [IN] YAC instance identifier

  call yac_ccleanup_instance_c ( yac_instance_id )

end subroutine yac_fcleanup_instance

! -------------------------- final -------------------------------------

subroutine yac_ffinalize ( )

  use mo_yac_finterface, dummy => yac_ffinalize

  implicit none

  interface
     subroutine yac_cfinalize_c () bind ( c, name='yac_cfinalize' )
     end subroutine yac_cfinalize_c
  end interface

  call yac_cfinalize_c ( )

end subroutine yac_ffinalize

subroutine yac_ffinalize_instance ( yac_instance_id )

  use mo_yac_finterface, dummy => yac_ffinalize_instance

  implicit none

  interface
     subroutine yac_cfinalize_instance_c ( yac_instance_id ) &
        bind ( c, name='yac_cfinalize_instance' )

        use, intrinsic :: iso_c_binding, only : c_int

        integer (kind=c_int), value :: yac_instance_id

     end subroutine yac_cfinalize_instance_c
  end interface

  integer, intent(in) :: yac_instance_id !< [IN] YAC instance identifier

  call yac_cfinalize_instance_c ( yac_instance_id )

end subroutine yac_ffinalize_instance

! -------------------------- version ----------------------------------

function yac_fget_version () result (version_string)

  use, intrinsic :: iso_c_binding, only : c_ptr
  use mo_yac_finterface, dummy => yac_fget_version
  use mo_yac_iso_c_helpers

  implicit none

  interface
     function yac_cget_version_c () bind ( c, name='yac_cget_version' )

     use, intrinsic :: iso_c_binding, only : c_ptr
     type(c_ptr) :: yac_cget_version_c

     end function yac_cget_version_c
  end interface

  type (c_ptr)                   :: c_string_ptr
  character (len=:), ALLOCATABLE :: version_string

  c_string_ptr = yac_cget_version_c()
  version_string = yac_internal_cptr2char(c_string_ptr)

end function yac_fget_version

! -------------------------- dates ------------------------------------

subroutine yac_fdef_datetime ( start_datetime, end_datetime )

  use, intrinsic :: iso_c_binding, only : c_null_char
  use mo_yac_finterface, dummy => yac_fdef_datetime
  use mo_yac_iso_c_helpers

  implicit none

  interface

     subroutine yac_cdef_datetime_c ( start_datetime, end_datetime ) &
    &      bind ( c, name='yac_cdef_datetime' )

       use, intrinsic :: iso_c_binding, only : c_char

       character ( kind=c_char), dimension(*) :: start_datetime
       character ( kind=c_char), dimension(*) :: end_datetime

     end subroutine yac_cdef_datetime_c

  end interface

  character(len=*), intent(in), optional :: start_datetime !< [IN] start datetime of job
  character(len=*), intent(in), optional :: end_datetime   !< [IN] end datetime of job

  integer :: index

  index = 0


  if (present(start_datetime)) then
    call yac_check_strlength ( start_datetime )
    index = index + 1
  end if

  if (present(end_datetime)) then
    call yac_check_strlength ( end_datetime )
    index = index + 2
  end if

  select case ( index )

    case ( 3 )
      call yac_cdef_datetime_c ( TRIM(start_datetime) // c_null_char, &
                                 TRIM(end_datetime)   // c_null_char )
    case ( 2 )
      call yac_cdef_datetime_c ( c_null_char, &
                                 TRIM(end_datetime)   // c_null_char )
    case ( 1 )
      call yac_cdef_datetime_c ( TRIM(start_datetime) // c_null_char, &
                                 c_null_char )
    end select

end subroutine yac_fdef_datetime

subroutine yac_fdef_datetime_instance ( &
  yac_instance_id, start_datetime, end_datetime )

  use, intrinsic :: iso_c_binding, only : c_null_char
  use mo_yac_finterface, dummy => yac_fdef_datetime_instance
  use mo_yac_iso_c_helpers

  implicit none

  interface

     subroutine yac_cdef_datetime_instance_c ( yac_instance_id, &
                                               start_datetime,  &
                                               end_datetime )   &
       bind ( c, name='yac_cdef_datetime_instance' )

       use, intrinsic :: iso_c_binding, only : c_char, c_int

       integer (kind=c_int), value            :: yac_instance_id
       character ( kind=c_char), dimension(*) :: start_datetime
       character ( kind=c_char), dimension(*) :: end_datetime

     end subroutine yac_cdef_datetime_instance_c

  end interface

  integer, intent(in)                    :: yac_instance_id !< [IN] YAC instance identifier
  character(len=*), intent(in), optional :: start_datetime  !< [IN] start datetime of job
  character(len=*), intent(in), optional :: end_datetime    !< [IN] end datetime of job

  integer :: index

  index = 0


  if (present(start_datetime)) then
    call yac_check_strlength ( start_datetime )
    index = index + 1
  end if

  if (present(end_datetime)) then
    call yac_check_strlength ( end_datetime )
    index = index + 2
  end if

  select case ( index )

    case ( 3 )
      call yac_cdef_datetime_instance_c ( yac_instance_id, &
                                          TRIM(start_datetime) // c_null_char, &
                                          TRIM(end_datetime)   // c_null_char )
    case ( 2 )
      call yac_cdef_datetime_instance_c ( yac_instance_id, &
                                          c_null_char, &
                                          TRIM(end_datetime)   // c_null_char )
    case ( 1 )
      call yac_cdef_datetime_instance_c ( yac_instance_id, &
                                          TRIM(start_datetime) // c_null_char, &
                                           c_null_char )
    end select

end subroutine yac_fdef_datetime_instance

subroutine yac_fdef_calendar ( calendar )

  use mo_yac_finterface, dummy => yac_fdef_calendar

  implicit none

  interface

      subroutine yac_cdef_calendar_c ( calendar ) &
        bind ( c, name='yac_cdef_calendar' )
        use, intrinsic :: iso_c_binding, only : c_int

        integer ( kind=c_int ), value :: calendar

      end subroutine yac_cdef_calendar_c

   end interface

  integer, intent(in) :: calendar

  call yac_cdef_calendar_c ( calendar )

end subroutine yac_fdef_calendar

! ------------------------ predef_comp ------------------------------------

SUBROUTINE yac_fpredef_comp ( comp_name, comp_id )

  use, intrinsic :: iso_c_binding, only : c_null_char
  use mo_yac_finterface, dummy => yac_fpredef_comp
  use mo_yac_iso_c_helpers

  implicit none

  INTERFACE

     SUBROUTINE yac_cpredef_comp_c ( comp_name, comp_id ) &
          bind ( c, name='yac_cpredef_comp' )

       use, intrinsic :: iso_c_binding, only : c_int, c_char

       character ( kind=c_char), dimension(*) :: comp_name
       integer   ( kind=c_int )               :: comp_id

     END SUBROUTINE yac_cpredef_comp_c

  END INTERFACE

  character(len=*), intent(in) :: comp_name !< [IN]  component name
  integer, intent(out)         :: comp_id   !< [OUT] returned handle to the component

  call yac_check_strlength ( comp_name )

  call yac_cpredef_comp_c ( TRIM(comp_name) // c_null_char, comp_id )

END SUBROUTINE yac_fpredef_comp

SUBROUTINE yac_fpredef_comp_instance ( yac_instance_id, comp_name, comp_id )

  use, intrinsic :: iso_c_binding, only : c_null_char
  use mo_yac_finterface, dummy => yac_fpredef_comp_instance
  use mo_yac_iso_c_helpers

  implicit none

  INTERFACE

     SUBROUTINE yac_cpredef_comp_instance_c ( yac_instance_id, &
          comp_name,       &
          comp_id )        &
          bind ( c, name='yac_cdef_comp_instance' )

       use, intrinsic :: iso_c_binding, only : c_int, c_char

       integer (kind=c_int), value            :: yac_instance_id
       character ( kind=c_char), dimension(*) :: comp_name
       integer   ( kind=c_int )               :: comp_id

     END SUBROUTINE yac_cpredef_comp_instance_c

  END INTERFACE

  integer, intent(in)          :: yac_instance_id !< [IN]  YAC instance identifier
  character(len=*), intent(in) :: comp_name       !< [IN]  component name
  integer, intent(out)         :: comp_id         !< [OUT] returned handle to the component

  call yac_check_strlength ( comp_name )

  call yac_cpredef_comp_instance_c( yac_instance_id,                &
       TRIM(comp_name) // c_null_char, &
       comp_id )

END SUBROUTINE yac_fpredef_comp_instance

! ------------------------ def_comp ------------------------------------

subroutine yac_fdef_comp ( comp_name, comp_id )

  use, intrinsic :: iso_c_binding, only : c_null_char
  use mo_yac_finterface, dummy => yac_fdef_comp
  use mo_yac_iso_c_helpers

  implicit none

  interface

     subroutine yac_cdef_comp_c ( comp_name, comp_id ) &
       bind ( c, name='yac_cdef_comp' )

       use, intrinsic :: iso_c_binding, only : c_int, c_char

       character ( kind=c_char), dimension(*) :: comp_name
       integer   ( kind=c_int )               :: comp_id

     end subroutine yac_cdef_comp_c

  end interface

  character(len=*), intent(in) :: comp_name !< [IN]  component name
  integer, intent(out)         :: comp_id   !< [OUT] returned handle to the component

  call yac_check_strlength ( comp_name )

  call yac_cdef_comp_c ( TRIM(comp_name) // c_null_char, comp_id )

end subroutine yac_fdef_comp

subroutine yac_fdef_comp_instance ( yac_instance_id, comp_name, comp_id )

  use, intrinsic :: iso_c_binding, only : c_null_char
  use mo_yac_finterface, dummy => yac_fdef_comp_instance
  use mo_yac_iso_c_helpers

  implicit none

  interface

     subroutine yac_cdef_comp_instance_c ( yac_instance_id, &
                                           comp_name,       &
                                           comp_id )        &
       bind ( c, name='yac_cdef_comp_instance' )

       use, intrinsic :: iso_c_binding, only : c_int, c_char

       integer (kind=c_int), value            :: yac_instance_id
       character ( kind=c_char), dimension(*) :: comp_name
       integer   ( kind=c_int )               :: comp_id

     end subroutine yac_cdef_comp_instance_c

  end interface

  integer, intent(in)          :: yac_instance_id !< [IN]  YAC instance identifier
  character(len=*), intent(in) :: comp_name       !< [IN]  component name
  integer, intent(out)         :: comp_id         !< [OUT] returned handle to the component

  call yac_check_strlength ( comp_name )

  call yac_cdef_comp_instance_c( yac_instance_id,                &
                                 TRIM(comp_name) // c_null_char, &
                                 comp_id )

end subroutine yac_fdef_comp_instance

! ------------------------ def_comps ------------------------------------

subroutine yac_fdef_comps ( comp_names, num_comps, comp_ids )

  use, intrinsic :: iso_c_binding, only : c_null_char, c_ptr, c_loc, c_char
  use mo_yac_finterface, dummy => yac_fdef_comps
  use mo_yac_iso_c_helpers

  implicit none

  interface

     subroutine yac_cdef_comps_c ( comp_names, num_comps, comp_ids ) &
           bind ( c, name='yac_cdef_comps' )

       use, intrinsic :: iso_c_binding, only : c_int, c_ptr

       type ( c_ptr )                :: comp_names(*)
       integer ( kind=c_int ), value :: num_comps
       integer ( kind=c_int )        :: comp_ids(*)

     end subroutine yac_cdef_comps_c

  end interface

  integer, intent(in)   :: num_comps !< [IN]  number of components
  character(kind=c_char, len=*), intent(in) :: &
                           comp_names(num_comps) !< [IN]  component names
  integer, intent(out)  :: comp_ids(num_comps) !< [OUT] returned handle to the components

  integer :: i, j
  character(kind=c_char), target :: comp_names_cpy(YAC_MAX_CHARLEN+1, num_comps)
  type(c_ptr) :: comp_name_ptrs(num_comps)

  comp_names_cpy = c_null_char

  do i = 1, num_comps
     call yac_check_strlength(comp_names(i))
    do j = 1, len_trim(comp_names(i))
      comp_names_cpy(j,i) = comp_names(i)(j:j)
    end do
    comp_name_ptrs(i) = c_loc(comp_names_cpy(1,i))
  end do

  call yac_cdef_comps_c ( comp_name_ptrs, num_comps, comp_ids )

end subroutine yac_fdef_comps

subroutine yac_fdef_comps_instance ( yac_instance_id, &
                                    comp_names,       &
                                    num_comps,        &
                                    comp_ids )

  use, intrinsic :: iso_c_binding, only : c_null_char, c_ptr, c_loc, c_char
  use mo_yac_finterface, dummy => yac_fdef_comps_instance
  use mo_yac_iso_c_helpers

  implicit none

  interface

     subroutine yac_cdef_comps_instance_c ( yac_instance_id, &
                                            comp_names,      &
                                            num_comps,       &
                                            comp_ids )       &
       bind ( c, name='yac_cdef_comps_instance' )

       use, intrinsic :: iso_c_binding, only : c_int, c_ptr

       integer (kind=c_int), value   :: yac_instance_id
       type ( c_ptr )                :: comp_names(*)
       integer ( kind=c_int ), value :: num_comps
       integer ( kind=c_int )        :: comp_ids(*)

     end subroutine yac_cdef_comps_instance_c

  end interface

  integer, intent(in)   :: yac_instance_id       !< [IN]  YAC instance identifier
  integer, intent(in)   :: num_comps             !< [IN]  number of components
  character(kind=c_char, len=*), intent(in) :: &
                           comp_names(num_comps) !< [IN]  component names
  integer, intent(out)  :: comp_ids(num_comps)   !< [OUT] returned handle to the components

  integer :: i, j
  character(kind=c_char), target :: comp_names_cpy(YAC_MAX_CHARLEN+1, num_comps)
  type(c_ptr) :: comp_name_ptrs(num_comps)

  comp_names_cpy = c_null_char

  do i = 1, num_comps
     call yac_check_strlength(comp_names(i))
    do j = 1, len_trim(comp_names(i))
      comp_names_cpy(j,i) = comp_names(i)(j:j)
    end do
    comp_name_ptrs(i) = c_loc(comp_names_cpy(1,i))
  end do

  call yac_cdef_comps_instance_c ( yac_instance_id, &
                                   comp_name_ptrs,  &
                                   num_comps,       &
                                   comp_ids )

end subroutine yac_fdef_comps_instance


! ------------------------- def_points ----------------------------------

subroutine yac_fdef_points_reg2d_real ( grid_id,       &
                                        nbr_points,    &
                                        location,      &
                                        x_points_real, &
                                        y_points_real, &
                                        point_id )

  use mo_yac_finterface, dummy => yac_fdef_points_reg2d_real

  implicit none

  integer, intent(in)  :: grid_id                      !< [IN]  grid identifier
  integer, intent(in)  :: nbr_points(2)                !< [IN]  number of points in x and y
  integer, intent(in)  :: location                     !< [IN]  location, one of center/edge/vertex
  real, intent(in)     :: x_points_real(nbr_points(1)) !< [IN]  longitudes of points
  real, intent(in)     :: y_points_real(nbr_points(2)) !< [IN]  latitudes of points
  integer, intent(out) :: point_id                     !< [OUT] point identifier

  double precision    ::  x_points(nbr_points(1))
  double precision    ::  y_points(nbr_points(2))

  x_points(:) = dble(x_points_real(:))
  y_points(:) = dble(y_points_real(:))

  call yac_fdef_points_reg2d_dble ( grid_id,    &
                                    nbr_points, &
                                    location,   &
                                    x_points,   &
                                    y_points,   &
                                    point_id )

end subroutine yac_fdef_points_reg2d_real

subroutine yac_fdef_points_reg2d_dble ( grid_id,    &
                                        nbr_points, &
                                        location,   &
                                        x_points,   &
                                        y_points,   &
                                        point_id )

  use mo_yac_finterface, dummy => yac_fdef_points_reg2d_dble

  implicit none

  interface

     subroutine yac_cdef_points_reg2d_c ( grid_id,      &
                                          nbr_points,   &
                                          location,     &
                                          x_points,     &
                                          y_points,     &
                                          point_id )    &
    &      bind ( c, name='yac_cdef_points_reg2d' )

       use, intrinsic :: iso_c_binding, only : c_int, c_double

       integer ( kind=c_int ), value :: grid_id
       integer ( kind=c_int )        :: nbr_points(2)
       integer ( kind=c_int ), value :: location

       real    ( kind=c_double )     :: x_points(nbr_points(1))
       real    ( kind=c_double )     :: y_points(nbr_points(2))

       integer ( kind=c_int )        :: point_id

     end subroutine yac_cdef_points_reg2d_c

  end interface

  integer, intent(in)  :: grid_id                         !< [IN]  grid identifier
  integer, intent(in)  :: nbr_points(2)                   !< [IN]  number of points in x and y
  integer, intent(in)  :: location                        !< [IN]  location, one of center/edge/vertex

  double precision, intent(in) :: x_points(nbr_points(1)) !< [IN]  longitudes of points
  double precision, intent(in) :: y_points(nbr_points(2)) !< [IN]  latitudes of points

  integer, intent(out) :: point_id                        !< [OUT] point identifier

  call yac_cdef_points_reg2d_c ( grid_id,    &
                                 nbr_points, &
                                 location,   &
                                 x_points,   &
                                 y_points,   &
                                 point_id )

end subroutine yac_fdef_points_reg2d_dble

subroutine yac_fdef_points_curve2d_real ( grid_id,       &
                                          nbr_points,    &
                                          location,      &
                                          x_points_real, &
                                          y_points_real, &
                                          point_id )

  use mo_yac_finterface, dummy => yac_fdef_points_curve2d_real

  implicit none

  integer, intent(in)  :: grid_id                      !< [IN]  grid identifier
  integer, intent(in)  :: nbr_points(2)                !< [IN]  number of points in x and y
  integer, intent(in)  :: location                     !< [IN]  location, one of center/edge/vertex
  real, intent(in)     :: &
    x_points_real(nbr_points(1),nbr_points(2))         !< [IN]  longitudes of points
  real, intent(in)     :: &
    y_points_real(nbr_points(1),nbr_points(2))         !< [IN]  latitudes of points
  integer, intent(out) :: point_id                     !< [OUT] point identifier

  double precision    ::  x_points(nbr_points(1),nbr_points(2))
  double precision    ::  y_points(nbr_points(1),nbr_points(2))

  x_points(:,:) = dble(x_points_real(:,:))
  y_points(:,:) = dble(y_points_real(:,:))

  call yac_fdef_points_curve2d_dble ( grid_id,    &
                                      nbr_points, &
                                      location,   &
                                      x_points,   &
                                      y_points,   &
                                      point_id )

end subroutine yac_fdef_points_curve2d_real

subroutine yac_fdef_points_curve2d_dble ( grid_id,    &
                                          nbr_points, &
                                          location,   &
                                          x_points,   &
                                          y_points,   &
                                          point_id )

  use mo_yac_finterface, dummy => yac_fdef_points_curve2d_dble

  implicit none

  interface

     subroutine yac_cdef_points_curve2d_c ( grid_id,    &
                                            nbr_points, &
                                            location,   &
                                            x_points,   &
                                            y_points,   &
                                            point_id )  &
           bind ( c, name='yac_cdef_points_curve2d' )

       use, intrinsic :: iso_c_binding, only : c_int, c_double

       integer ( kind=c_int ), value :: grid_id
       integer ( kind=c_int )        :: nbr_points(2)
       integer ( kind=c_int ), value :: location

       real    ( kind=c_double )     :: x_points(nbr_points(1),nbr_points(2))
       real    ( kind=c_double )     :: y_points(nbr_points(1),nbr_points(2))

       integer ( kind=c_int )        :: point_id

     end subroutine yac_cdef_points_curve2d_c

  end interface

  integer, intent(in)  :: grid_id                         !< [IN]  grid identifier
  integer, intent(in)  :: nbr_points(2)                   !< [IN]  number of points in x and y
  integer, intent(in)  :: location                        !< [IN]  location, one of center/edge/vertex

  double precision, intent(in) :: &
    x_points(nbr_points(1),nbr_points(2))                 !< [IN]  longitudes of points
  double precision, intent(in) :: &
    y_points(nbr_points(1),nbr_points(2))                 !< [IN]  latitudes of points

  integer, intent(out) :: point_id                        !< [OUT] point identifier

  call yac_cdef_points_curve2d_c ( grid_id,    &
                                   nbr_points, &
                                   location,   &
                                   x_points,   &
                                   y_points,   &
                                   point_id )

end subroutine yac_fdef_points_curve2d_dble

subroutine yac_fdef_points_unstruct_real ( grid_id,       &
                                           nbr_points,    &
                                           location,      &
                                           x_points_real, &
                                           y_points_real, &
                                           point_id )

  use mo_yac_finterface, dummy => yac_fdef_points_unstruct_real

  implicit none

  integer, intent(in)  :: grid_id                   !< [IN]  grid identifier
  integer, intent(in)  :: nbr_points                !< [IN]  number of points
  integer, intent(in)  :: location                  !< [IN]  location, one of center/edge/vertex

  real, intent(in)     :: x_points_real(nbr_points) !< [IN]  longitudes of points
  real, intent(in)     :: y_points_real(nbr_points) !< [IN]  latitudes of points

  integer, intent(out) :: point_id                  !< [OUT] point identifier

  double precision    ::  x_points(nbr_points)
  double precision    ::  y_points(nbr_points)

  x_points(:) = dble(x_points_real(:))
  y_points(:) = dble(y_points_real(:))

  call yac_fdef_points_unstruct_dble ( grid_id,    &
                                       nbr_points, &
                                       location,   &
                                       x_points,   &
                                       y_points,   &
                                       point_id )

end subroutine yac_fdef_points_unstruct_real

subroutine yac_fdef_points_unstruct_dble ( grid_id,    &
                                           nbr_points, &
                                           location,   &
                                           x_points,   &
                                           y_points,   &
                                           point_id )

  use mo_yac_finterface, dummy => yac_fdef_points_unstruct_dble

  implicit none

  interface

     subroutine yac_cdef_points_unstruct_c ( grid_id,      &
                                             nbr_points,   &
                                             location,     &
                                             x_points,     &
                                             y_points,     &
                                             point_id )    &
           bind ( c, name='yac_cdef_points_unstruct' )

       use, intrinsic :: iso_c_binding, only : c_int, c_double

       integer (kind=c_int), value :: grid_id
       integer (kind=c_int), value :: nbr_points
       integer (kind=c_int), value :: location

       real    (kind=c_double)     :: x_points(nbr_points)
       real    (kind=c_double)     :: y_points(nbr_points)

       integer (kind=c_int)        :: point_id

     end subroutine yac_cdef_points_unstruct_c

  end interface

  integer, intent(in)  :: grid_id                      !< [IN]  grid identifier
  integer, intent(in)  :: nbr_points                   !< [IN]  number of points
  integer, intent(in)  :: location                     !< [IN]  location, one of center/edge/vertex

  double precision, intent(in) :: x_points(nbr_points) !< [IN]  longitudes of points
  double precision, intent(in) :: y_points(nbr_points) !< [IN]  latitudes of points

  integer, intent(out) :: point_id                     !< [OUT] grid identifier

  call yac_cdef_points_unstruct_c  ( grid_id,    &
                                     nbr_points, &
                                     location,   &
                                     x_points,   &
                                     y_points,   &
                                     point_id )

end subroutine yac_fdef_points_unstruct_dble

! ------------------------- def_grid -------------------------------

!> Definition of a non-uniform unstructured grid (cells have
!! varying numbers of vertices)
!! @param[in]  grid_name             grid name
!! @param[in]  nbr_vertices          number of vertices
!! @param[in]  nbr_cells             number of cells
!! @param[in]  nbr_connections       total size of cell_to_vertex
!! @param[in]  nbr_vertices_per_cell number of vertices for each cell
!! @param[in]  x_vertices_real       longitudes of vertices
!! @param[in]  y_vertices_real       latitudes of vertices
!! @param[in]  cell_to_vertex_in     connectivity between vertices and cells\n
!!                                   (the vertex indices per cell have to be
!!                                   in clockwise or counterclockwise ordering)
!! @param[out] grid_id               grid identifier
!! @param[in]  use_ll_edges          use lonlat edges
!!
!! @remark If (use_ll_edges == .TRUE.) YAC will check all edges of the grid
!!         and determine whether they are on circles of longitudes
!!         (same x coordinate) or latitudes (same y coordinate). An edge that
!!         does not fullfill this condition will cause an error.
subroutine yac_fdef_grid_nonuniform_real ( grid_name,             &
                                           nbr_vertices,          &
                                           nbr_cells,             &
                                           nbr_connections,       &
                                           nbr_vertices_per_cell, &
                                           x_vertices_real,       &
                                           y_vertices_real,       &
                                           cell_to_vertex_in,     &
                                           grid_id,               &
                                           use_ll_edges)

  use mo_yac_finterface, dummy => yac_fdef_grid_nonuniform_real

  implicit none

  character(len=*), intent(in) :: grid_name
  integer, intent(in)  :: nbr_vertices
  integer, intent(in)  :: nbr_cells
  integer, intent(in)  :: nbr_connections
  integer, intent(in)  :: nbr_vertices_per_cell(nbr_cells)

  real, intent(in)     :: x_vertices_real(nbr_vertices)
  real, intent(in)     :: y_vertices_real(nbr_vertices)

  integer, intent(in)  :: cell_to_vertex_in(nbr_connections)

  integer, intent(out) :: grid_id

  logical, optional, intent(in) :: use_ll_edges

  double precision     :: x_vertices(nbr_vertices)
  double precision     :: y_vertices(nbr_vertices)

  x_vertices(:) = dble(x_vertices_real(:))
  y_vertices(:) = dble(y_vertices_real(:))

  call yac_fdef_grid_nonuniform_dble ( grid_name,             &
                                       nbr_vertices,          &
                                       nbr_cells,             &
                                       nbr_connections,       &
                                       nbr_vertices_per_cell, &
                                       x_vertices,            &
                                       y_vertices,            &
                                       cell_to_vertex_in,     &
                                       grid_id,               &
                                       use_ll_edges )

end subroutine yac_fdef_grid_nonuniform_real

!> Definition of a non-uniform unstructured grid (cells have
!! varying numbers of vertices)
!! @param[in]  grid_name             grid name
!! @param[in]  nbr_vertices          number of vertices
!! @param[in]  nbr_cells             number of cells
!! @param[in]  nbr_connections       total size of cell_to_vertex
!! @param[in]  nbr_vertices_per_cell number of vertices for each cell
!! @param[in]  x_vertices            longitudes of vertices
!! @param[in]  y_vertices            latitudes of vertices
!! @param[in]  cell_to_vertex_in     connectivity between vertices and cells\n
!!                                   (the vertex indices per cell have to be in
!!                                   clockwise or counterclockwise ordering)
!! @param[out] grid_id               grid identifier
!! @param[in]  use_ll_edges          use lonlat edges
!!
!! @remark If (use_ll_edges == .TRUE.) YAC will check all edges of the grid
!!         and determine whether they are on circles of longitudes
!!         (same x coordinate) or latitudes (same y coordinate). An edge that
!!         does not fullfill this condition will cause an error.
subroutine yac_fdef_grid_nonuniform_dble ( grid_name,             &
                                           nbr_vertices,          &
                                           nbr_cells,             &
                                           nbr_connections,       &
                                           nbr_vertices_per_cell, &
                                           x_vertices,            &
                                           y_vertices,            &
                                           cell_to_vertex_in,     &
                                           grid_id,               &
                                           use_ll_edges )

  use, intrinsic :: iso_c_binding, only : c_null_char
  use mo_yac_finterface, dummy => yac_fdef_grid_nonuniform_dble

  implicit none

  interface

     subroutine yac_cdef_grid_unstruct_c ( grid_name,             &
                                           nbr_vertices,          &
                                           nbr_cells,             &
                                           nbr_vertices_per_cell, &
                                           x_vertices,            &
                                           y_vertices,            &
                                           cell_to_vertex,        &
                                           grid_id )              &
           bind ( c, name='yac_cdef_grid_unstruct' )

       use, intrinsic :: iso_c_binding, only : c_int, c_char, c_double

       character ( kind=c_char ), dimension(*) :: grid_name
       integer ( kind=c_int    ), value        :: nbr_vertices
       integer ( kind=c_int    ), value        :: nbr_cells
       integer ( kind=c_int    )               :: nbr_vertices_per_cell(nbr_cells)
       real    ( kind=c_double )               :: x_vertices(nbr_vertices)
       real    ( kind=c_double )               :: y_vertices(nbr_vertices)
       integer ( kind=c_int    )               :: cell_to_vertex(nbr_cells,nbr_vertices)
       integer ( kind=c_int    )               :: grid_id

     end subroutine yac_cdef_grid_unstruct_c

     subroutine yac_cdef_grid_unstruct_ll_c ( grid_name,             &
                                              nbr_vertices,          &
                                              nbr_cells,             &
                                              nbr_vertices_per_cell, &
                                              x_vertices,            &
                                              y_vertices,            &
                                              cell_to_vertex,        &
                                              grid_id )              &
           bind ( c, name='yac_cdef_grid_unstruct_ll' )

       use, intrinsic :: iso_c_binding, only : c_int, c_char, c_double

       character ( kind=c_char ), dimension(*) :: grid_name
       integer ( kind=c_int    ), value        :: nbr_vertices
       integer ( kind=c_int    ), value        :: nbr_cells
       integer ( kind=c_int    )               :: nbr_vertices_per_cell(nbr_cells)
       real    ( kind=c_double )               :: x_vertices(nbr_vertices)
       real    ( kind=c_double )               :: y_vertices(nbr_vertices)
       integer ( kind=c_int    )               :: cell_to_vertex(nbr_cells,nbr_vertices)
       integer ( kind=c_int    )               :: grid_id

     end subroutine yac_cdef_grid_unstruct_ll_c

  end interface

  character(len=*), intent(in) :: grid_name
  integer, intent(in)          :: nbr_vertices
  integer, intent(in)          :: nbr_cells
  integer, intent(in)          :: nbr_connections
  integer, intent(in)          :: nbr_vertices_per_cell(nbr_cells)

  double precision, intent(in) :: x_vertices(nbr_vertices)
  double precision, intent(in) :: y_vertices(nbr_vertices)

  integer, intent(in)          :: cell_to_vertex_in(nbr_connections)

  integer, intent(out)         :: grid_id

  logical, optional, intent(in) :: use_ll_edges


  integer                      :: cell_to_vertex(nbr_connections)

  logical :: use_ll_edges_

  YAC_FASSERT(ALL(cell_to_vertex_in > 0), "ERROR(yac_fdef_grid_nonuniform_dble): all entries of cell_to_vertex have to be > 0")

  YAC_FASSERT(ALL(cell_to_vertex_in <= nbr_vertices), "ERROR(yac_fdef_grid_nonuniform_dble): all entries of cell_to_vertex have to be <= nbr_vertices")

  cell_to_vertex(:) = cell_to_vertex_in(:) - 1

  if (present(use_ll_edges)) then
    use_ll_edges_ = use_ll_edges
  else
    use_ll_edges_ = .false.
  end if

  if (use_ll_edges_) then
    call yac_cdef_grid_unstruct_ll_c  ( TRIM(grid_name) // c_null_char, &
                                        nbr_vertices,                   &
                                        nbr_cells,                      &
                                        nbr_vertices_per_cell,          &
                                        x_vertices,                     &
                                        y_vertices,                     &
                                        cell_to_vertex,                 &
                                        grid_id )
  else
    call yac_cdef_grid_unstruct_c  ( TRIM(grid_name) // c_null_char, &
                                     nbr_vertices,                   &
                                     nbr_cells,                      &
                                     nbr_vertices_per_cell,          &
                                     x_vertices,                     &
                                     y_vertices,                     &
                                     cell_to_vertex,                 &
                                     grid_id )
  end if

end subroutine yac_fdef_grid_nonuniform_dble

!> Definition of a uniform unstructured grid (all cells have the
!! number of vertices)
!! @param[in]  grid_name                grid name
!! @param[in]  nbr_vertices             number of vertices
!! @param[in]  nbr_cells                number of cells
!! @param[in]  nbr_vertices_per_cell_in number of vertices for each cell
!! @param[in]  x_vertices_real          longitudes of vertices
!! @param[in]  y_vertices_real          latitudes of vertices
!! @param[in]  cell_to_vertex_in        connectivity between vertices and cells\n
!!                                      (the vertex indices per cell have to be
!!                                      in clockwise or counterclockwise ordering)
!! @param[out] grid_id                  grid identifier
!! @param[in]  use_ll_edges             use lonlat edges
!!
!! @remark If (use_ll_edges == .TRUE.) YAC will check all edges of the grid
!!         and determine whether they are on circles of longitudes
!!         (same x coordinate) or latitudes (same y coordinate). An edge that
!!         does not fullfill this condition will cause an error.
subroutine yac_fdef_grid_unstruct_real ( grid_name,                &
                                         nbr_vertices,             &
                                         nbr_cells,                &
                                         nbr_vertices_per_cell_in, &
                                         x_vertices_real,          &
                                         y_vertices_real,          &
                                         cell_to_vertex_in,        &
                                         grid_id,                  &
                                         use_ll_edges )

  use mo_yac_finterface, dummy => yac_fdef_grid_unstruct_real

  implicit none

  character(len=*), intent(in) :: grid_name
  integer, intent(in)  :: nbr_vertices
  integer, intent(in)  :: nbr_cells
  integer, intent(in)  :: nbr_vertices_per_cell_in

  real, intent(in)     :: x_vertices_real(nbr_vertices)
  real, intent(in)     :: y_vertices_real(nbr_vertices)

  integer, intent(in)  :: cell_to_vertex_in(nbr_vertices_per_cell_in,nbr_cells)

  integer, intent(out)   :: grid_id

  logical, optional, intent(in) :: use_ll_edges

  double precision       :: x_vertices(nbr_vertices)
  double precision       :: y_vertices(nbr_vertices)

  x_vertices(:) = dble(x_vertices_real(:))
  y_vertices(:) = dble(y_vertices_real(:))

  call yac_fdef_grid_unstruct_dble ( grid_name,                &
                                     nbr_vertices,             &
                                     nbr_cells,                &
                                     nbr_vertices_per_cell_in, &
                                     x_vertices,               &
                                     y_vertices,               &
                                     cell_to_vertex_in,        &
                                     grid_id,                  &
                                     use_ll_edges )

end subroutine yac_fdef_grid_unstruct_real

!> Definition of a uniform unstructured grid (all cells have the
!! number of vertices)
!! @param[in]  grid_name                grid name
!! @param[in]  nbr_vertices             number of vertices
!! @param[in]  nbr_cells                number of cells
!! @param[in]  nbr_vertices_per_cell_in number of vertices for each cell
!! @param[in]  x_vertices               longitudes of vertices
!! @param[in]  y_vertices               latitudes of vertices
!! @param[in]  cell_to_vertex_in        connectivity between vertices and cells\n
!!                                      (the vertex indices per cell have to be
!!                                      in clockwise or counterclockwise ordering)
!! @param[out] grid_id                  grid identifier
!! @param[in]  use_ll_edges             use lonlat edges
!!
!! @remark If (use_ll_edges == .TRUE.) YAC will check all edges of the grid
!!         and determine whether they are on circles of longitudes
!!         (same x coordinate) or latitudes (same y coordinate). An edge that
!!         does not fullfill this condition will cause an error.
subroutine yac_fdef_grid_unstruct_dble ( grid_name,                &
                                         nbr_vertices,             &
                                         nbr_cells,                &
                                         nbr_vertices_per_cell_in, &
                                         x_vertices,               &
                                         y_vertices,               &
                                         cell_to_vertex_in,        &
                                         grid_id,                  &
                                         use_ll_edges )

  use, intrinsic :: iso_c_binding, only : c_null_char
  use mo_yac_finterface, dummy => yac_fdef_grid_unstruct_dble

  implicit none

  interface

     subroutine yac_cdef_grid_unstruct_c ( grid_name,             &
                                           nbr_vertices,          &
                                           nbr_cells,             &
                                           nbr_vertices_per_cell, &
                                           x_vertices,            &
                                           y_vertices,            &
                                           cell_to_vertex,        &
                                           grid_id )              &
           bind ( c, name='yac_cdef_grid_unstruct' )

       use, intrinsic :: iso_c_binding, only : c_int, c_char, c_double

       character ( kind=c_char ), dimension(*) :: grid_name
       integer ( kind=c_int    ), value        :: nbr_vertices
       integer ( kind=c_int    ), value        :: nbr_cells
       integer ( kind=c_int    )               :: nbr_vertices_per_cell(nbr_cells)
       real    ( kind=c_double )               :: x_vertices(nbr_vertices)
       real    ( kind=c_double )               :: y_vertices(nbr_vertices)
       integer ( kind=c_int    )               :: cell_to_vertex(nbr_cells,nbr_vertices)
       integer ( kind=c_int    )               :: grid_id

     end subroutine yac_cdef_grid_unstruct_c

     subroutine yac_cdef_grid_unstruct_ll_c ( grid_name,             &
                                              nbr_vertices,          &
                                              nbr_cells,             &
                                              nbr_vertices_per_cell, &
                                              x_vertices,            &
                                              y_vertices,            &
                                              cell_to_vertex,        &
                                              grid_id )              &
           bind ( c, name='yac_cdef_grid_unstruct_ll' )

       use, intrinsic :: iso_c_binding, only : c_int, c_char, c_double

       character ( kind=c_char ), dimension(*) :: grid_name
       integer ( kind=c_int    ), value        :: nbr_vertices
       integer ( kind=c_int    ), value        :: nbr_cells
       integer ( kind=c_int    )               :: nbr_vertices_per_cell(nbr_cells)
       real    ( kind=c_double )               :: x_vertices(nbr_vertices)
       real    ( kind=c_double )               :: y_vertices(nbr_vertices)
       integer ( kind=c_int    )               :: cell_to_vertex(nbr_cells,nbr_vertices)
       integer ( kind=c_int    )               :: grid_id

     end subroutine yac_cdef_grid_unstruct_ll_c

  end interface

  character(len=*), intent(in) :: grid_name
  integer, intent(in)          :: nbr_vertices
  integer, intent(in)          :: nbr_cells
  integer, intent(in)          :: nbr_vertices_per_cell_in

  double precision, intent(in) :: x_vertices(nbr_vertices)
  double precision, intent(in) :: y_vertices(nbr_vertices)

  integer, intent(in)  :: cell_to_vertex_in( &
                            nbr_vertices_per_cell_in, &
                            nbr_cells)

  integer, intent(out)         :: grid_id

  logical, optional, intent(in) :: use_ll_edges

  integer                      :: nbr_vertices_per_cell(nbr_cells)
  integer                      :: cell_to_vertex(nbr_vertices_per_cell_in,nbr_cells)
  logical                      :: use_ll_edges_

  nbr_vertices_per_cell(:) = nbr_vertices_per_cell_in

  YAC_FASSERT(ALL(cell_to_vertex_in > 0), "ERROR(yac_fdef_grid_unstruct_dble): all entries of cell_to_vertex have to be > 0")

  YAC_FASSERT(ALL(cell_to_vertex_in <= nbr_vertices), "ERROR(yac_fdef_grid_unstruct_dble): all entries of cell_to_vertex have to be <= nbr_vertices")

  cell_to_vertex(:,:) = cell_to_vertex_in(:,:) - 1

  if (present(use_ll_edges)) then
    use_ll_edges_ = use_ll_edges
  else
    use_ll_edges_ = .false.
  end if

  if (use_ll_edges_) then
    call yac_cdef_grid_unstruct_ll_c  ( TRIM(grid_name) // c_null_char, &
                                        nbr_vertices,                   &
                                        nbr_cells,                      &
                                        nbr_vertices_per_cell,          &
                                        x_vertices,                     &
                                        y_vertices,                     &
                                        cell_to_vertex,                 &
                                        grid_id )
  else
    call yac_cdef_grid_unstruct_c  ( TRIM(grid_name) // c_null_char, &
                                     nbr_vertices,                   &
                                     nbr_cells,                      &
                                     nbr_vertices_per_cell,          &
                                     x_vertices,                     &
                                     y_vertices,                     &
                                     cell_to_vertex,                 &
                                     grid_id )
  end if

end subroutine yac_fdef_grid_unstruct_dble

!> Definition of a 2d curvilinear grid
!! @param[in]  grid_name       grid name
!! @param[in]  nbr_vertices    number of cells in each dimension
!! @param[in]  cyclic          cyclic behavior of cells in each dimension
!! @param[in]  x_vertices_real longitudes of vertices
!! @param[in]  y_vertices_real latitudes of vertices
!! @param[out] grid_id         grid identifier
subroutine yac_fdef_grid_curve2d_real ( grid_name,       &
                                        nbr_vertices,    &
                                        cyclic,          &
                                        x_vertices_real, &
                                        y_vertices_real, &
                                        grid_id )

  use mo_yac_finterface, dummy => yac_fdef_grid_curve2d_real

  implicit none

  character(len=*), intent(in) :: grid_name
  integer, intent(in)  :: nbr_vertices(2)
  integer, intent(in)  :: cyclic(2)
  real, intent(in)     :: &
    x_vertices_real(nbr_vertices(1),nbr_vertices(2))
  real, intent(in)     :: &
    y_vertices_real(nbr_vertices(1),nbr_vertices(2))
  integer, intent(out) :: grid_id

  double precision     ::  x_vertices(nbr_vertices(1),nbr_vertices(2))
  double precision     ::  y_vertices(nbr_vertices(1),nbr_vertices(2))

  x_vertices(:,:) = dble(x_vertices_real(:,:))
  y_vertices(:,:) = dble(y_vertices_real(:,:))

  call yac_fdef_grid_curve2d_dble ( grid_name,    &
                                    nbr_vertices, &
                                    cyclic,       &
                                    x_vertices,   &
                                    y_vertices,   &
                                    grid_id )

end subroutine yac_fdef_grid_curve2d_real

!> Definition of a 2d curvilinear grid
!! @param[in]  grid_name    grid name
!! @param[in]  nbr_vertices number of cells in each dimension
!! @param[in]  cyclic       cyclic behavior of cells in each dimension
!! @param[in]  x_vertices   longitudes of vertices
!! @param[in]  y_vertices   latitudes of vertices
!! @param[out] grid_id      grid identifier
subroutine yac_fdef_grid_curve2d_dble ( grid_name,    &
                                        nbr_vertices, &
                                        cyclic,       &
                                        x_vertices,   &
                                        y_vertices,   &
                                        grid_id )

  use, intrinsic :: iso_c_binding, only : c_null_char
  use mo_yac_finterface, dummy => yac_fdef_grid_curve2d_dble

  implicit none

  interface

     subroutine yac_cdef_grid_curve2d_c ( grid_name,    &
                                          nbr_vertices, &
                                          cyclic,       &
                                          x_vertices,   &
                                          y_vertices,   &
                                          grid_id )     &
           bind ( c, name='yac_cdef_grid_curve2d' )

       use, intrinsic :: iso_c_binding, only : c_int, c_char, c_double

       character ( kind=c_char ), dimension(*) :: grid_name
       integer ( kind=c_int    )               :: nbr_vertices(2)
       integer ( kind=c_int    )               :: cyclic(2)
       real    ( kind=c_double )               :: x_vertices(nbr_vertices(1),nbr_vertices(2))
       real    ( kind=c_double )               :: y_vertices(nbr_vertices(1),nbr_vertices(2))
       integer ( kind=c_int    )               :: grid_id

     end subroutine yac_cdef_grid_curve2d_c

  end interface

  character(len=*), intent(in) :: grid_name
  integer, intent(in)          :: nbr_vertices(2)
  integer, intent(in)          :: cyclic(2)
  double precision, intent(in) :: &
    x_vertices(nbr_vertices(1),nbr_vertices(2))
  double precision, intent(in) :: &
    y_vertices(nbr_vertices(1),nbr_vertices(2))
  integer, intent(out)         :: grid_id

  call yac_cdef_grid_curve2d_c ( TRIM(grid_name) // c_null_char, &
                                 nbr_vertices,                   &
                                 cyclic,                         &
                                 x_vertices,                     &
                                 y_vertices,                     &
                                 grid_id )

end subroutine yac_fdef_grid_curve2d_dble

!> Definition of a 2d regular grid
!! @param[in]  grid_name       grid name
!! @param[in]  nbr_vertices    number of cells in each dimension
!! @param[in]  cyclic          cyclic behavior of cells in each dimension
!! @param[in]  x_vertices_real longitudes of vertices
!! @param[in]  y_vertices_real latitudes of vertices
!! @param[out] grid_id         grid identifier
subroutine yac_fdef_grid_reg2d_real ( grid_name,       &
                                      nbr_vertices,    &
                                      cyclic,          &
                                      x_vertices_real, &
                                      y_vertices_real, &
                                      grid_id )

  use mo_yac_finterface, dummy => yac_fdef_grid_reg2d_real

  implicit none

  character(len=*), intent(in) :: grid_name
  integer, intent(in)  :: nbr_vertices(2)
  integer, intent(in)  :: cyclic(2)
  real, intent(in)     :: x_vertices_real(nbr_vertices(1))
  real, intent(in)     :: y_vertices_real(nbr_vertices(2))
  integer, intent(out) :: grid_id

  double precision     ::  x_vertices(nbr_vertices(1))
  double precision     ::  y_vertices(nbr_vertices(2))

  x_vertices(:) = dble(x_vertices_real(:))
  y_vertices(:) = dble(y_vertices_real(:))

  call yac_fdef_grid_reg2d_dble ( grid_name,    &
                                  nbr_vertices, &
                                  cyclic,       &
                                  x_vertices,   &
                                  y_vertices,   &
                                  grid_id )

end subroutine yac_fdef_grid_reg2d_real

!> Definition of a 2d regular grid
!! @param[in]  grid_name    grid name
!! @param[in]  nbr_vertices number of cells in each dimension
!! @param[in]  cyclic       cyclic behavior of cells in each dimension
!! @param[in]  x_vertices   longitudes of vertices
!! @param[in]  y_vertices   latitudes of vertices
!! @param[out] grid_id      grid identifier
subroutine yac_fdef_grid_reg2d_dble ( grid_name,    &
                                      nbr_vertices, &
                                      cyclic,       &
                                      x_vertices,   &
                                      y_vertices,   &
                                      grid_id )

  use, intrinsic :: iso_c_binding, only : c_null_char
  use mo_yac_finterface, dummy => yac_fdef_grid_reg2d_dble

  implicit none

  interface

     subroutine yac_cdef_grid_reg2d_c ( grid_name,    &
                                        nbr_vertices, &
                                        cyclic,       &
                                        x_vertices,   &
                                        y_vertices,   &
                                        grid_id )     &
           bind ( c, name='yac_cdef_grid_reg2d' )

       use, intrinsic :: iso_c_binding, only : c_int, c_char, c_double

       character ( kind=c_char ), dimension(*) :: grid_name
       integer ( kind=c_int    )               :: nbr_vertices(2)
       integer ( kind=c_int    )               :: cyclic(2)
       real    ( kind=c_double )               :: x_vertices(nbr_vertices(1))
       real    ( kind=c_double )               :: y_vertices(nbr_vertices(2))
       integer ( kind=c_int    )               :: grid_id

     end subroutine yac_cdef_grid_reg2d_c

  end interface

  character(len=*), intent(in) :: grid_name
  integer, intent(in)          :: nbr_vertices(2)
  integer, intent(in)          :: cyclic(2)
  double precision, intent(in) :: x_vertices(nbr_vertices(1))
  double precision, intent(in) :: y_vertices(nbr_vertices(2))
  integer, intent(out)         :: grid_id

  call yac_cdef_grid_reg2d_c ( TRIM(grid_name) // c_null_char, &
                               nbr_vertices,                   &
                               cyclic,                         &
                               x_vertices,                     &
                               y_vertices,                     &
                               grid_id )

end subroutine yac_fdef_grid_reg2d_dble

! ------------------------- set global_index ------------------------------

subroutine yac_fset_global_index ( global_index, &
                                   location,     &
                                   grid_id )

  use mo_yac_finterface, dummy => yac_fset_global_index

  implicit none

  interface

     subroutine yac_cset_global_index_c ( global_index, &
                                          location,     &
                                          grid_id )     &
           bind ( c, name='yac_cset_global_index' )

       use, intrinsic :: iso_c_binding, only : c_int

       integer ( kind=c_int )        :: global_index(*)
       integer ( kind=c_int ), value :: location
       integer ( kind=c_int ), value :: grid_id

     end subroutine yac_cset_global_index_c

  end interface

  integer, intent(in)  :: global_index(*) !< [IN] global index
  integer, intent(in)  :: location        !< [IN] location
  integer, intent(in)  :: grid_id         !< [IN] point set indentifier

  call yac_cset_global_index_c ( global_index, &
                                 location,     &
                                 grid_id )

end subroutine yac_fset_global_index

! ------------------------- set core_mask ------------------------------

subroutine yac_fset_core_lmask ( is_core,  &
                                 location, &
                                 grid_id )

  use mo_yac_finterface, dummy => yac_fset_core_lmask
  use, intrinsic :: iso_c_binding, only : c_size_t

  implicit none

  logical, intent(in)  :: is_core(*) !< [IN] core flag
                                     !< false, cell/vertex/edge is halo
                                     !< true, cell/vertex/edge is core
  integer, intent(in)  :: location   !< [IN] location
  integer, intent(in)  :: grid_id    !< [IN] point set indentifier

  integer (kind=c_size_t) :: i, count
  integer, allocatable :: int_is_core(:)

  interface

     function yac_get_grid_size_c ( location, &
                                    grid_id ) &
           result ( grid_size )               &
           bind ( c, name='yac_get_grid_size' )

       use, intrinsic :: iso_c_binding, only : c_int, c_size_t

       integer ( kind=c_int ), value :: location
       integer ( kind=c_int ), value :: grid_id
       integer ( kind=c_size_t) :: grid_size

     end function yac_get_grid_size_c

  end interface

  count = yac_get_grid_size_c(location, grid_id)
  allocate(int_is_core(count))

  do i = 1, count
     if ( is_core(i) ) then
        int_is_core(i) = 1
     else
        int_is_core(i) = 0
     endif
  enddo

  call yac_fset_core_imask ( int_is_core, &
                             location,    &
                             grid_id )

end subroutine yac_fset_core_lmask

subroutine yac_fset_core_imask ( is_core,  &
                                 location, &
                                 grid_id )

  use mo_yac_finterface, dummy => yac_fset_core_imask

  implicit none

  interface

     subroutine yac_cset_core_mask_c ( mask,        &
                                       location,    &
                                       grid_id )    &
           bind ( c, name='yac_cset_core_mask' )

       use, intrinsic :: iso_c_binding, only : c_int

       integer ( kind=c_int )        :: mask(*)
       integer ( kind=c_int ), value :: location
       integer ( kind=c_int ), value :: grid_id

     end subroutine yac_cset_core_mask_c

  end interface


  integer, intent(in)  :: is_core(*) !< [IN] core flag
                                     !< 0, cell/vertex/edge is halo
                                     !< 1, cell/vertex/edge is core
  integer, intent(in)  :: location   !< [IN] location
  integer, intent(in)  :: grid_id    !< [IN] point set indentifier

  call yac_cset_core_mask_c ( is_core,  &
                              location, &
                              grid_id )

end subroutine yac_fset_core_imask

! ------------------------- set_mask ------------------------------

subroutine yac_fset_lmask ( is_valid,     &
                            points_id )

  use mo_yac_finterface, dummy => yac_fset_lmask
  use, intrinsic :: iso_c_binding, only : c_size_t

  implicit none

  logical, intent(in)  :: is_valid(*)  !< [IN] logical mask
                                       !< false, point is masked out
                                       !< true, point is valid
  integer, intent(in)  :: points_id    !< [IN] point set indentifier

  integer ( kind=c_size_t) :: i, count
  integer, allocatable :: int_is_valid(:)

  interface

     function yac_get_points_size_c ( points_id ) &
           result ( points_size )                 &
           bind ( c, name='yac_get_points_size' )

       use, intrinsic :: iso_c_binding, only : c_int, c_size_t

       integer ( kind=c_int ), value :: points_id
       integer ( kind=c_size_t) :: points_size

     end function yac_get_points_size_c

  end interface

  count = yac_get_points_size_c(points_id)
  allocate(int_is_valid(count))

  do i = 1, count
    int_is_valid(i) = MERGE(1,0,is_valid(i))
  enddo

  call yac_fset_imask ( int_is_valid, points_id )

end subroutine yac_fset_lmask

subroutine yac_fset_imask ( is_valid, &
                            points_id )

  use mo_yac_finterface, dummy => yac_fset_imask

  implicit none

  interface

     subroutine yac_cset_mask_c ( is_valid,   &
                                  points_id ) &
           bind ( c, name='yac_cset_mask' )

       use, intrinsic :: iso_c_binding, only : c_int

       integer ( kind=c_int )        :: is_valid(*)
       integer ( kind=c_int ), value :: points_id

     end subroutine yac_cset_mask_c

  end interface

  integer, intent(in)  :: is_valid(*)  !< [IN] integer mask
                                       !< 0, point is masked out
                                       !< 1, point is valid
  integer, intent(in)  :: points_id    !< [IN] point set indentifier

  call yac_cset_mask_c ( is_valid, points_id )

end subroutine yac_fset_imask

! ------------------------- def_mask ------------------------------

subroutine yac_fdef_lmask ( grid_id,    &
                            nbr_points, &
                            location,   &
                            is_valid,   &
                            mask_id )

  use mo_yac_finterface, dummy => yac_fdef_lmask

  implicit none

  integer, intent(in)  :: grid_id     !< [IN] grid identifier
  integer, intent(in)  :: nbr_points  !< [IN] number of points
  integer, intent(in)  :: location    !< [IN] location, one of center/edge/vertex
  logical, intent(in)  :: is_valid(*) !< [IN] logical mask
                                      !< false, point is masked out
                                      !< true, point is valid
  integer, intent(out) :: mask_id     !< [OUT] mask identifier

  integer :: i
  integer, allocatable :: int_is_valid(:)

  allocate(int_is_valid(nbr_points))

  do i = 1, nbr_points
     int_is_valid(i) = MERGE(1,0,is_valid(i))
  enddo

  call yac_fdef_imask ( grid_id,      &
                        nbr_points,   &
                        location,     &
                        int_is_valid, &
                        mask_id )

end subroutine yac_fdef_lmask

subroutine yac_fdef_imask ( grid_id,    &
                            nbr_points, &
                            location,   &
                            is_valid,   &
                            mask_id )

  use mo_yac_finterface, dummy => yac_fdef_imask

  implicit none

  interface

     subroutine yac_cdef_mask_c ( grid_id,    &
                                  nbr_points, &
                                  location,   &
                                  is_valid,   &
                                  mask_id )   &
           bind ( c, name='yac_cdef_mask' )

       use, intrinsic :: iso_c_binding, only : c_int

       integer ( kind=c_int ), value :: grid_id
       integer ( kind=c_int ), value :: nbr_points
       integer ( kind=c_int ), value :: location
       integer ( kind=c_int )        :: is_valid(*)
       integer ( kind=c_int )        :: mask_id

     end subroutine yac_cdef_mask_c

  end interface

  integer, intent(in)  :: grid_id     !< [IN] grid identifier
  integer, intent(in)  :: nbr_points  !< [IN] number of points
  integer, intent(in)  :: location    !< [IN] location, one of center/edge/vertex
  integer, intent(in)  :: is_valid(*) !< [IN] logical mask
                                      !< false, point is masked out
                                      !< true, point is valid
  integer, intent(out) :: mask_id     !< [OUT] mask identifier

  call yac_cdef_mask_c ( grid_id,    &
                         nbr_points, &
                         location,   &
                         is_valid,   &
                         mask_id )

end subroutine yac_fdef_imask

! ------------------------- def_mask_named ------------------------

subroutine yac_fdef_lmask_named ( grid_id,    &
                                  nbr_points, &
                                  location,   &
                                  is_valid,   &
                                  name,       &
                                  mask_id )

  use mo_yac_finterface, dummy => yac_fdef_lmask_named

  implicit none

  integer, intent(in)  :: grid_id      !< [IN] grid identifier
  integer, intent(in)  :: nbr_points   !< [IN] number of points
  integer, intent(in)  :: location     !< [IN] location, one of center/edge/vertex
  logical, intent(in)  :: is_valid(*)  !< [IN] logical mask
                                       !< false, point is masked out
                                       !< true, point is valid
  character(len=*), intent(in) :: name !< [IN] name of the mask
  integer, intent(out) :: mask_id      !< [OUT] mask identifier

  integer :: i
  integer, allocatable :: int_is_valid(:)

  allocate(int_is_valid(nbr_points))

  do i = 1, nbr_points
     int_is_valid(i) = MERGE(1,0,is_valid(i))
  enddo

  call yac_fdef_imask_named ( grid_id,      &
                              nbr_points,   &
                              location,     &
                              int_is_valid, &
                              name,         &
                              mask_id )

end subroutine yac_fdef_lmask_named

subroutine yac_fdef_imask_named ( grid_id,    &
                                  nbr_points, &
                                  location,   &
                                  is_valid,   &
                                  name,       &
                                  mask_id )

  use, intrinsic :: iso_c_binding, only : c_null_char
  use mo_yac_finterface, dummy => yac_fdef_imask_named
  use mo_yac_iso_c_helpers

  implicit none

  interface

     subroutine yac_cdef_mask_named_c ( grid_id,    &
                                        nbr_points, &
                                        location,   &
                                        is_valid,   &
                                        name,       &
                                        mask_id )   &
           bind ( c, name='yac_cdef_mask_named' )

       use, intrinsic :: iso_c_binding, only : c_int, c_char

       integer ( kind=c_int ), value           :: grid_id
       integer ( kind=c_int ), value           :: nbr_points
       integer ( kind=c_int ), value           :: location
       integer ( kind=c_int )                  :: is_valid(*)
       character ( kind=c_char ), dimension(*) :: name
       integer ( kind=c_int )                  :: mask_id

     end subroutine yac_cdef_mask_named_c

  end interface

  integer, intent(in)  :: grid_id      !< [IN] grid identifier
  integer, intent(in)  :: nbr_points   !< [IN] number of points
  integer, intent(in)  :: location     !< [IN] location, one of center/edge/vertex
  integer, intent(in)  :: is_valid(*)  !< [IN] logical mask
                                       !< false, point is masked out
                                       !< true, point is valid
  character(len=*), intent(in) :: name !< [IN] name of the mask
  integer, intent(out) :: mask_id      !< [OUT] mask identifier

  call yac_check_strlength ( name )

  call yac_cdef_mask_named_c ( grid_id,                   &
                               nbr_points,                &
                               location,                  &
                               is_valid,                  &
                               TRIM(name) // c_null_char, &
                               mask_id )

end subroutine yac_fdef_imask_named

! ----------------------------- def_field -------------------------------

subroutine yac_fdef_field ( field_name,      &
                            component_id,    &
                            point_ids,       &
                            num_pointsets,   &
                            collection_size, &
                            timestep,        &
                            time_unit,       &
                            field_id )

  use, intrinsic :: iso_c_binding, only : c_null_char
  use mo_yac_finterface, dummy => yac_fdef_field
  use mo_yac_iso_c_helpers

  implicit none

  interface

     subroutine yac_cdef_field_c ( field_name,      &
                                   component_id,    &
                                   point_ids,       &
                                   num_pointsets,   &
                                   collection_size, &
                                   timestep,        &
                                   time_unit,       &
                                   field_id )       &
           bind ( c, name='yac_cdef_field' )

       use, intrinsic :: iso_c_binding, only : c_int, c_char

       character ( kind=c_char ), dimension(*) :: field_name !< [IN] short name of the field
       integer ( kind=c_int ), value :: component_id         !< [IN] component identifier
       integer ( kind=c_int )        :: point_ids(*)         !< [IN] point identifier
       integer ( kind=c_int ), value :: num_pointsets        !< [IN] number of pointsets per grid
       integer ( kind=c_int ), value :: collection_size      !< [IN] collection_size
       character ( kind=c_char ), dimension(*) :: timestep   !< [IN] timestep
       integer ( kind=c_int ), value :: time_unit            !< [IN] unit of timestep
       integer ( kind=c_int )        :: field_id             !< [OUT] returned field handle

     end subroutine yac_cdef_field_c

  end interface

  !
  ! Definition of coupling fields
  !
  character(len=*), intent (in) :: field_name      !< [IN] short name of the field
  integer, intent (in)          :: component_id    !< [IN] component identifier
  integer, intent (in)          :: point_ids(*)    !< [IN] point identifier
  integer, intent (in)          :: num_pointsets   !< [IN] number of pointsets per grid
  integer, intent (in)          :: collection_size !< [IN] collection size
  character(len=*), intent (in) :: timestep        !< [IN] timestep
  integer, intent (in)          :: time_unit       !< [IN] unit of timestep
  integer, intent (out)         :: field_id        !< [OUT] returned field handle

  call yac_check_strlength ( field_name )
  call yac_check_strlength ( timestep )


  call yac_cdef_field_c ( TRIM(field_name) // c_null_char, &
                          component_id,                    &
                          point_ids,                       &
                          num_pointsets,                   &
                          collection_size,                 &
                          TRIM(timestep) // c_null_char,   &
                          time_unit,                       &
                          field_id )

end subroutine yac_fdef_field

! ----------------------------- def_field_mask---------------------------

subroutine yac_fdef_field_mask ( field_name,      &
                                 component_id,    &
                                 point_ids,       &
                                 mask_ids,        &
                                 num_pointsets,   &
                                 collection_size, &
                                 timestep,        &
                                 time_unit,       &
                                 field_id )

  use, intrinsic :: iso_c_binding, only : c_null_char
  use mo_yac_finterface, dummy => yac_fdef_field_mask
  use mo_yac_iso_c_helpers

  implicit none

  interface

     subroutine yac_cdef_field_mask_c ( field_name,      &
                                        component_id,    &
                                        point_ids,       &
                                        mask_ids,        &
                                        num_pointsets,   &
                                        collection_size, &
                                        timestep,        &
                                        time_unit,       &
                                        field_id )       &
           bind ( c, name='yac_cdef_field_mask' )

       use, intrinsic :: iso_c_binding, only : c_int, c_char

       character ( kind=c_char ), dimension(*) :: field_name
       integer ( kind=c_int ), value :: component_id
       integer ( kind=c_int )        :: point_ids(*)
       integer ( kind=c_int )        :: mask_ids(*)
       integer ( kind=c_int ), value :: num_pointsets
       integer ( kind=c_int ), value :: collection_size
       character ( kind=c_char ), dimension(*) :: timestep
       integer ( kind=c_int ), value :: time_unit
       integer ( kind=c_int )        :: field_id

     end subroutine yac_cdef_field_mask_c

  end interface


  !
  ! Definition of coupling fields
  !
  character(len=*), intent (in) :: field_name      !< [IN] short name of the field
  integer, intent (in)          :: component_id    !< [IN] component identifier
  integer, intent (in)          :: point_ids(*)    !< [IN] point identifier
  integer, intent (in)          :: mask_ids(*)     !< [IN] mask identifier
  integer, intent (in)          :: num_pointsets   !< [IN] number of pointsets per grid
  integer, intent (in)          :: collection_size !< [IN] collection size
  character(len=*), intent (in) :: timestep        !< [IN] timestep
  integer, intent (in)          :: time_unit       !< [IN] unit of timestep
  integer, intent (out)         :: field_id        !< [OUT] returned field handle

  call yac_check_strlength ( field_name )
  call yac_check_strlength ( timestep )

  call yac_cdef_field_mask_c ( TRIM(field_name) // c_null_char, &
                               component_id,                    &
                               point_ids,                       &
                               mask_ids,                        &
                               num_pointsets,                   &
                               collection_size,                 &
                               TRIM(timestep) // c_null_char,   &
                               time_unit,                       &
                               field_id )

end subroutine yac_fdef_field_mask

! -----------------------------------------------------------------------

subroutine yac_fcheck_field_dimensions( field_id,          &
                                        collection_size,   &
                                        num_interp_fields, &
                                        interp_field_sizes )

  use mo_yac_finterface, dummy => yac_fcheck_field_dimensions
  use, intrinsic :: iso_c_binding, only : c_int

  implicit none

  interface

     subroutine yac_ccheck_field_dimensions_c ( field_id,            &
                                                collection_size,     &
                                                num_interp_fields,   &
                                                interp_field_sizes ) &
       bind ( c, name='yac_ccheck_field_dimensions' )

       use, intrinsic :: iso_c_binding, only : c_int

       integer ( kind=c_int ), value :: field_id
       integer ( kind=c_int ), value :: collection_size
       integer ( kind=c_int ), value :: num_interp_fields
       integer ( kind=c_int ), dimension(*) :: interp_field_sizes

     end subroutine yac_ccheck_field_dimensions_c

  end interface

  integer, intent (in) :: field_id                              !<[IN] field handle
  integer, intent (in) :: collection_size                       !<[IN] collection size
  integer, intent (in) :: num_interp_fields                     !<[IN] number of interpolation fields
                                                                !!     (number of pointsets)
  integer, intent (in) :: interp_field_sizes(num_interp_fields) !<[IN] data size of each
                                                                !!     interpolation field

  call yac_ccheck_field_dimensions_c(field_id,          &
                                     collection_size,   &
                                     num_interp_fields, &
                                     interp_field_sizes)

end subroutine yac_fcheck_field_dimensions

! ---------------------------------- put --------------------------------
!>
!! @param[in]  field_id        field identifier
!! @param[in]  nbr_hor_points  number of horizontal points
!! @param[in]  nbr_pointsets   number of point sets
!! @param[in]  collection_size number of vertical level or bundles
!! @param[in]  send_field      returned field
!! @param[out] info            returned info
!! @param[out] ierror          returned error
subroutine yac_fput_real ( field_id,         &
                           nbr_hor_points,   &
                           nbr_pointsets,    &
                           collection_size,  &
                           send_field,       &
                           info,             &
                           ierror )

  use mo_yac_finterface, dummy => yac_fput_real
  use mo_yac_real_to_dble_utils

  implicit none

  integer, intent (in)  :: field_id
  integer, intent (in)  :: nbr_hor_points
  integer, intent (in)  :: nbr_pointsets
  integer, intent (in)  :: collection_size
  real, intent (in)     :: send_field(nbr_hor_points, &
                                      nbr_pointsets,  &
                                      collection_size)
  integer, intent (out) :: info
  integer, intent (out) :: ierror

  double precision      :: send_field_dble(nbr_hor_points, &
                                           nbr_pointsets,  &
                                           collection_size)
  integer :: i

  call yac_fcheck_field_dimensions(           &
    field_id, collection_size, nbr_pointsets, &
    (/(nbr_hor_points,i=1,nbr_pointsets)/) )

  call send_field_to_dble(field_id,        &
                          nbr_hor_points,  &
                          nbr_pointsets,   &
                          collection_size, &
                          send_field,      &
                          send_field_dble)

  call yac_fput_dble ( field_id,        &
                       nbr_hor_points,  &
                       nbr_pointsets,   &
                       collection_size, &
                       send_field_dble, &
                       info,            &
                       ierror )

end subroutine yac_fput_real

!>
!! @param[in]  field_id        field identifier
!! @param[in]  nbr_hor_points  number of horizontal points
!! @param[in]  nbr_pointsets   number of point sets
!! @param[in]  collection_size number of vertical level or bundles
!! @param[in]  send_field      send field
!! @param[in]  send_frac_mask  fractional mask
!! @param[out] info            returned info
!! @param[out] ierror          returned error
subroutine yac_fput_frac_real ( field_id,         &
                                nbr_hor_points,   &
                                nbr_pointsets,    &
                                collection_size,  &
                                send_field,       &
                                send_frac_mask,   &
                                info,             &
                                ierror )

  use mo_yac_finterface, dummy => yac_fput_frac_real
  use mo_yac_real_to_dble_utils

  implicit none

  integer, intent (in)  :: field_id
  integer, intent (in)  :: nbr_hor_points
  integer, intent (in)  :: nbr_pointsets
  integer, intent (in)  :: collection_size
  real, intent (in)     :: send_field(nbr_hor_points, &
                                      nbr_pointsets,  &
                                      collection_size)
  real, intent (in)     :: send_frac_mask(nbr_hor_points, &
                                          nbr_pointsets,  &
                                          collection_size)
  integer, intent (out) :: info
  integer, intent (out) :: ierror

  double precision      :: send_field_dble(nbr_hor_points, &
                                           nbr_pointsets,  &
                                           collection_size)
  double precision      :: send_frac_mask_dble(nbr_hor_points, &
                                               nbr_pointsets,  &
                                               collection_size)
  integer :: i

  call yac_fcheck_field_dimensions(           &
    field_id, collection_size, nbr_pointsets, &
    (/(nbr_hor_points,i=1,nbr_pointsets)/) )

  call send_field_to_dble(field_id,        &
                          nbr_hor_points,  &
                          nbr_pointsets,   &
                          collection_size, &
                          send_field,      &
                          send_field_dble, &
                          send_frac_mask,  &
                          send_frac_mask_dble)

  call yac_fput_frac_dble ( field_id,            &
                            nbr_hor_points,      &
                            nbr_pointsets,       &
                            collection_size,     &
                            send_field_dble,     &
                            send_frac_mask_dble, &
                            info,                &
                            ierror )

end subroutine yac_fput_frac_real

!>
!! @param[in]  field_id        field identifier
!! @param[in]  nbr_pointsets   number of point sets
!! @param[in]  collection_size number of vertical level or bundles
!! @param[in]  send_field      send field
!! @param[out] info            returned info
!! @param[out] ierror          returned error
subroutine yac_fput_real_ptr ( field_id,        &
                               nbr_pointsets,   &
                               collection_size, &
                               send_field,      &
                               info,            &
                               ierror )

  use mo_yac_finterface, dummy => yac_fput_real_ptr
  use mo_yac_real_to_dble_utils

  implicit none

  integer, intent (in)            :: field_id
  integer, intent (in)            :: nbr_pointsets
  integer, intent (in)            :: collection_size
  type(yac_real_ptr), intent (in) :: send_field(nbr_pointsets, &
                                                collection_size)
  integer, intent (out)           :: info
  integer, intent (out)           :: ierror

  integer :: i, j
  type(yac_dble_ptr)              :: send_field_dble(nbr_pointsets, &
                                                     collection_size)

  call yac_fcheck_field_dimensions(           &
    field_id, collection_size, nbr_pointsets, &
    (/(SIZE(send_field(i, 1)%p),i=1,nbr_pointsets)/) )

  call send_field_to_dble_ptr(field_id,        &
                              nbr_pointsets,   &
                              collection_size, &
                              send_field,      &
                              send_field_dble)

  call yac_fput_dble_ptr ( field_id,        &
                           nbr_pointsets,   &
                           collection_size, &
                           send_field_dble, &
                           info,            &
                           ierror )

  do i = 1, collection_size
    do j = 1, nbr_pointsets
      deallocate(send_field_dble(j, i)%p)
    end do
  end do

end subroutine yac_fput_real_ptr

!>
!! @param[in]  field_id        field identifier
!! @param[in]  nbr_pointsets   number of point sets
!! @param[in]  collection_size number of vertical level or bundles
!! @param[in]  send_field      send field
!! @param[in]  send_frac_mask  fractional mask
!! @param[out] info            returned info
!! @param[out] ierror          returned error
subroutine yac_fput_frac_real_ptr ( field_id,        &
                                    nbr_pointsets,   &
                                    collection_size, &
                                    send_field,      &
                                    send_frac_mask,  &
                                    info,            &
                                    ierror )

  use mo_yac_finterface, dummy => yac_fput_frac_real_ptr
  use mo_yac_real_to_dble_utils

  implicit none

  integer, intent (in)            :: field_id
  integer, intent (in)            :: nbr_pointsets
  integer, intent (in)            :: collection_size
  type(yac_real_ptr), intent (in) :: send_field(nbr_pointsets, &
                                                collection_size)
  type(yac_real_ptr), intent (in) :: send_frac_mask(nbr_pointsets, &
                                                    collection_size)
  integer, intent (out)           :: info
  integer, intent (out)           :: ierror

  integer :: i, j
  type(yac_dble_ptr)              :: send_field_dble(nbr_pointsets, &
                                                     collection_size)
  type(yac_dble_ptr)              :: send_frac_mask_dble(nbr_pointsets, &
                                                         collection_size)

  call yac_fcheck_field_dimensions(           &
    field_id, collection_size, nbr_pointsets, &
    (/(SIZE(send_field(i, 1)%p),i=1,nbr_pointsets)/) )

  call send_field_to_dble_ptr(field_id,        &
                              nbr_pointsets,   &
                              collection_size, &
                              send_field,      &
                              send_field_dble, &
                              send_frac_mask,  &
                              send_frac_mask_dble)

  call yac_fput_frac_dble_ptr ( field_id,            &
                                nbr_pointsets,       &
                                collection_size,     &
                                send_field_dble,     &
                                send_frac_mask_dble, &
                                info,                &
                                ierror )

  do i = 1, collection_size
    do j = 1, nbr_pointsets
      deallocate(send_field_dble(j, i)%p)
      deallocate(send_frac_mask_dble(j, i)%p)
    end do
  end do

end subroutine yac_fput_frac_real_ptr

!>
!! @param[in]  field_id        field identifier
!! @param[in]  nbr_hor_points  number of horizontal points
!! @param[in]  collection_size number of vertical level or bundles
!! @param[in]  send_field      send field
!! @param[out] info            returned info
!! @param[out] ierror          returned error
subroutine yac_fput_single_pointset_real ( field_id,        &
                                           nbr_hor_points,  &
                                           collection_size, &
                                           send_field,      &
                                           info,            &
                                           ierror )

  use mo_yac_finterface, dummy => yac_fput_single_pointset_real
  use mo_yac_real_to_dble_utils

  implicit none

  integer, intent (in)          :: field_id
  integer, intent (in)          :: nbr_hor_points
  integer, intent (in)          :: collection_size
  real, intent (in)             :: send_field(nbr_hor_points, &
                                              collection_size)
  integer, intent (out)         :: info
  integer, intent (out)         :: ierror

  double precision              :: send_field_dble(nbr_hor_points, &
                                                   collection_size)

  call yac_fcheck_field_dimensions( &
    field_id, collection_size, 1, (/nbr_hor_points/) )

  call send_field_to_dble_single(field_id,        &
                                 nbr_hor_points,  &
                                 collection_size, &
                                 send_field,      &
                                 send_field_dble)

  call yac_fput_single_pointset_dble ( field_id,        &
                                       nbr_hor_points,  &
                                       collection_size, &
                                       send_field_dble, &
                                       info,            &
                                       ierror )

end subroutine yac_fput_single_pointset_real

!>
!! @param[in]  field_id        field identifier
!! @param[in]  nbr_hor_points  number of horizontal points
!! @param[in]  collection_size number of vertical level or bundles
!! @param[in]  send_field      send field
!! @param[in]  send_frac_mask  fractional mask
!! @param[out] info            returned info
!! @param[out] ierror          returned error
subroutine yac_fput_frac_single_pointset_real ( field_id,        &
                                                nbr_hor_points,  &
                                                collection_size, &
                                                send_field,      &
                                                send_frac_mask,  &
                                                info,            &
                                                ierror )

  use mo_yac_finterface, dummy => yac_fput_frac_single_pointset_real
  use mo_yac_real_to_dble_utils

  implicit none

  integer, intent (in)          :: field_id
  integer, intent (in)          :: nbr_hor_points
  integer, intent (in)          :: collection_size
  real, intent (in)             :: send_field(nbr_hor_points, &
                                              collection_size)
  real, intent (in)             :: send_frac_mask(nbr_hor_points, &
                                                  collection_size)
  integer, intent (out)         :: info
  integer, intent (out)         :: ierror

  double precision              :: send_field_dble(nbr_hor_points, &
                                                   collection_size)
  double precision              :: send_frac_mask_dble(nbr_hor_points, &
                                                       collection_size)
  call yac_fcheck_field_dimensions( &
    field_id, collection_size, 1, (/nbr_hor_points/) )

  call send_field_to_dble_single(field_id,        &
                                 nbr_hor_points,  &
                                 collection_size, &
                                 send_field,      &
                                 send_field_dble, &
                                 send_frac_mask,  &
                                 send_frac_mask_dble)

  call yac_fput_frac_single_pointset_dble ( field_id,            &
                                            nbr_hor_points,      &
                                            collection_size,     &
                                            send_field_dble,     &
                                            send_frac_mask_dble, &
                                            info,                &
                                            ierror )

end subroutine yac_fput_frac_single_pointset_real

!>
!! @param[in]  field_id        field identifier
!! @param[in]  nbr_hor_points  number of horizontal points
!! @param[in]  collection_size number of vertical level or bundles
!! @param[in]  nbr_pointsets   number of point sets
!! @param[in]  send_field      send field
!! @param[out] info            returned info
!! @param[out] ierror          returned error
subroutine yac_fput_dble ( field_id,        &
                           nbr_hor_points,  &
                           nbr_pointsets,   &
                           collection_size, &
                           send_field,      &
                           info,            &
                           ierror )

  use mo_yac_finterface, dummy => yac_fput_dble

  implicit none

  interface

     subroutine yac_cput__c ( field_id,         &
                              collection_size,  &
                              send_field,       &
                              info,             &
                              ierror ) bind ( c, name='yac_cput_' )

       use, intrinsic :: iso_c_binding, only : c_int, c_double

       integer ( kind=c_int ), value :: field_id
       integer ( kind=c_int ), value :: collection_size
       real    ( kind=c_double )     :: send_field(*)
       integer ( kind=c_int )        :: info
       integer ( kind=c_int )        :: ierror

     end subroutine yac_cput__c

  end interface

  integer, intent (in)          :: field_id
  integer, intent (in)          :: nbr_hor_points
  integer, intent (in)          :: nbr_pointsets
  integer, intent (in)          :: collection_size
  double precision, intent (in) :: send_field(nbr_hor_points, &
                                              nbr_pointsets,  &
                                              collection_size)
  integer, intent (out)         :: info
  integer, intent (out)         :: ierror

  integer :: i

  call yac_fcheck_field_dimensions(           &
    field_id, collection_size, nbr_pointsets, &
    (/(nbr_hor_points,i=1,nbr_pointsets)/) )
  YAC_FASSERT(is_contiguous(send_field), "ERROR(yac_fput_dble): send_field is not contiguous")

  call yac_cput__c ( field_id,         &
                     collection_size,  &
                     send_field,       &
                     info,             &
                     ierror )

end subroutine yac_fput_dble

!>
!! @param[in]  field_id        field identifier
!! @param[in]  nbr_hor_points  number of horizontal points
!! @param[in]  collection_size number of vertical level or bundles
!! @param[in]  nbr_pointsets   number of point sets
!! @param[in]  send_field      send field
!! @param[in]  send_frac_mask  fractional mask
!! @param[out] info            returned info
!! @param[out] ierror          returned error
subroutine yac_fput_frac_dble ( field_id,        &
                                nbr_hor_points,  &
                                nbr_pointsets,   &
                                collection_size, &
                                send_field,      &
                                send_frac_mask,  &
                                info,            &
                                ierror )

  use mo_yac_finterface, dummy => yac_fput_frac_dble

  implicit none

  interface

     subroutine yac_cput_frac__c ( field_id,         &
                                   collection_size,  &
                                   send_field,       &
                                   send_frac_mask,   &
                                   info,             &
                                   ierror ) &
       bind ( c, name='yac_cput_frac_' )

       use, intrinsic :: iso_c_binding, only : c_int, c_double

       integer ( kind=c_int ), value :: field_id
       integer ( kind=c_int ), value :: collection_size
       real    ( kind=c_double )     :: send_field(*)
       real    ( kind=c_double )     :: send_frac_mask(*)
       integer ( kind=c_int )        :: info
       integer ( kind=c_int )        :: ierror

     end subroutine yac_cput_frac__c

  end interface

  integer, intent (in)          :: field_id
  integer, intent (in)          :: nbr_hor_points
  integer, intent (in)          :: nbr_pointsets
  integer, intent (in)          :: collection_size
  double precision, intent (in) :: send_field(nbr_hor_points, &
                                              nbr_pointsets,  &
                                              collection_size)
  double precision, intent (in) :: send_frac_mask(nbr_hor_points, &
                                                  nbr_pointsets,  &
                                                  collection_size)
  integer, intent (out)         :: info
  integer, intent (out)         :: ierror

  integer :: i

  call yac_fcheck_field_dimensions(           &
    field_id, collection_size, nbr_pointsets, &
    (/(nbr_hor_points,i=1,nbr_pointsets)/) )
  YAC_FASSERT(is_contiguous(send_field), "ERROR(yac_fput_frac_dble): send_field is not contiguous")

  call yac_cput_frac__c ( field_id,         &
                          collection_size,  &
                          send_field,       &
                          send_frac_mask,   &
                          info,             &
                          ierror )

end subroutine yac_fput_frac_dble

!>
!! @param[in]  field_id        field identifier
!! @param[in]  nbr_pointsets   number of point sets
!! @param[in]  collection_size number of vertical level or bundles
!! @param[in]  send_field      send field
!! @param[out] info            returned info
!! @param[out] ierror          returned error
subroutine yac_fput_dble_ptr ( field_id,        &
                               nbr_pointsets,   &
                               collection_size, &
                               send_field,      &
                               info,            &
                               ierror )

  use mo_yac_finterface, dummy => yac_fput_dble_ptr
  use iso_c_binding, only: c_ptr, c_loc

  implicit none

  interface

     subroutine yac_cput_ptr__c ( field_id,        &
                                  collection_size, &
                                  send_field,      &
                                  info,            &
                                  ierror )         &
         bind ( c, name='yac_cput_ptr_' )

       use, intrinsic :: iso_c_binding, only : c_int, c_ptr

       integer ( kind=c_int ), value :: field_id
       integer ( kind=c_int ), value :: collection_size
       type(c_ptr)                   :: send_field(*)
       integer ( kind=c_int )        :: info
       integer ( kind=c_int )        :: ierror

     end subroutine yac_cput_ptr__c

  end interface

  integer, intent (in)     :: field_id
  integer, intent (in)     :: nbr_pointsets
  integer, intent (in)     :: collection_size
  type(yac_dble_ptr), intent (in) :: send_field(nbr_pointsets, collection_size)
  integer, intent (out)    :: info
  integer, intent (out)    :: ierror

  integer :: i, j
  type(c_ptr) :: send_field_(nbr_pointsets, collection_size)

  call yac_fcheck_field_dimensions(           &
    field_id, collection_size, nbr_pointsets, &
    (/(SIZE(send_field(i, 1)%p),i=1,nbr_pointsets)/) )

  do i = 1, collection_size
    do j = 1, nbr_pointsets
      YAC_FASSERT(is_contiguous(send_field(j, i)%p), "ERROR(yac_fput_dble_ptr): send_field is not contiguous")
      send_field_(j, i) = c_loc(send_field(j, i)%p(1))
    end do
  end do

  call yac_cput_ptr__c ( field_id,        &
                         collection_size, &
                         send_field_,     &
                         info,            &
                         ierror )

end subroutine yac_fput_dble_ptr

!>
!! @param[in]  field_id        field identifier
!! @param[in]  nbr_pointsets   number of point sets
!! @param[in]  collection_size number of vertical level or bundles
!! @param[in]  send_field      send field
!! @param[in]  send_frac_mask  fractional mask
!! @param[out] info            returned info
!! @param[out] ierror          returned error
subroutine yac_fput_frac_dble_ptr ( field_id,        &
                                    nbr_pointsets,   &
                                    collection_size, &
                                    send_field,      &
                                    send_frac_mask,  &
                                    info,            &
                                    ierror )

  use mo_yac_finterface, dummy => yac_fput_frac_dble_ptr
  use iso_c_binding, only: c_ptr, c_loc

  implicit none

  interface

     subroutine yac_cput_frac_ptr__c ( field_id,        &
                                       collection_size, &
                                       send_field,      &
                                       send_frac_mask,  &
                                       info,            &
                                       ierror )         &
         bind ( c, name='yac_cput_frac_ptr_' )

       use, intrinsic :: iso_c_binding, only : c_int, c_ptr

       integer ( kind=c_int ), value :: field_id
       integer ( kind=c_int ), value :: collection_size
       type(c_ptr)                   :: send_field(*)
       type(c_ptr)                   :: send_frac_mask(*)
       integer ( kind=c_int )        :: info
       integer ( kind=c_int )        :: ierror

     end subroutine yac_cput_frac_ptr__c

  end interface

  integer, intent (in)     :: field_id
  integer, intent (in)     :: nbr_pointsets
  integer, intent (in)     :: collection_size
  type(yac_dble_ptr), intent (in) :: send_field(nbr_pointsets, collection_size)
  type(yac_dble_ptr), intent (in) :: send_frac_mask(nbr_pointsets, collection_size)
  integer, intent (out)    :: info
  integer, intent (out)    :: ierror

  integer :: i, j
  type(c_ptr) :: send_field_(nbr_pointsets, collection_size)
  type(c_ptr) :: send_frac_mask_(nbr_pointsets, collection_size)

  call yac_fcheck_field_dimensions(           &
    field_id, collection_size, nbr_pointsets, &
    (/(SIZE(send_field(i, 1)%p),i=1,nbr_pointsets)/) )

  do i = 1, collection_size
    do j = 1, nbr_pointsets
      YAC_FASSERT(is_contiguous(send_field(j, i)%p), "ERROR(yac_fput_frac_dble_ptr): send_field is not contiguous")
      send_field_(j, i) = c_loc(send_field(j, i)%p(1))
      YAC_FASSERT(is_contiguous(send_frac_mask(j, i)%p), "ERROR(yac_fput_frac_dble_ptr): send_frac_mask is not contiguous")
      send_frac_mask_(j, i) = c_loc(send_frac_mask(j, i)%p(1))
    end do
  end do

  call yac_cput_frac_ptr__c ( field_id,        &
                              collection_size, &
                              send_field_,     &
                              send_frac_mask_, &
                              info,            &
                              ierror )

end subroutine yac_fput_frac_dble_ptr

!>
!! @param[in]  field_id        field identifier
!! @param[in]  nbr_hor_points  number of horizontal points
!! @param[in]  collection_size number of vertical level or bundles
!! @param[in]  send_field      send field
!! @param[out] info            returned info
!! @param[out] ierror          returned error
subroutine yac_fput_single_pointset_dble ( field_id,        &
                                           nbr_hor_points,  &
                                           collection_size, &
                                           send_field,      &
                                           info,            &
                                           ierror )

  use mo_yac_finterface, dummy => yac_fput_single_pointset_dble

  implicit none

  interface

     subroutine yac_cput__c ( field_id,        &
                              collection_size, &
                              send_field,      &
                              info,            &
                              ierror ) bind ( c, name='yac_cput_' )

       use, intrinsic :: iso_c_binding, only : c_int, c_double

       integer ( kind=c_int ), value :: field_id
       integer ( kind=c_int ), value :: collection_size
       real    ( kind=c_double )     :: send_field(*)
       integer ( kind=c_int )        :: info
       integer ( kind=c_int )        :: ierror

     end subroutine yac_cput__c

  end interface

  integer, intent (in)          :: field_id
  integer, intent (in)          :: nbr_hor_points
  integer, intent (in)          :: collection_size
  double precision, intent (in) :: send_field(nbr_hor_points, &
                                              collection_size)
  integer, intent (out)         :: info
  integer, intent (out)         :: ierror

  call yac_fcheck_field_dimensions( &
    field_id, collection_size, 1, (/nbr_hor_points/) )
  YAC_FASSERT(is_contiguous(send_field), "ERROR(yac_fput_single_pointset_dble): send_field is not contiguous")

  call yac_cput__c ( field_id,        &
                     collection_size, &
                     send_field,      &
                     info,            &
                     ierror )

end subroutine yac_fput_single_pointset_dble

!>
!! @param[in]  field_id        field identifier
!! @param[in]  nbr_hor_points  number of horizontal points
!! @param[in]  collection_size number of vertical level or bundles
!! @param[in]  send_field      send field
!! @param[in]  send_frac_mask  fractional mask
!! @param[out] info            returned info
!! @param[out] ierror          returned error
subroutine yac_fput_frac_single_pointset_dble ( field_id,        &
                                                nbr_hor_points,  &
                                                collection_size, &
                                                send_field,      &
                                                send_frac_mask,  &
                                                info,            &
                                                ierror )

  use mo_yac_finterface, dummy => yac_fput_frac_single_pointset_dble

  implicit none

  interface

     subroutine yac_cput_frac__c ( field_id,        &
                                   collection_size, &
                                   send_field,      &
                                   send_frac_mask,  &
                                   info,            &
                                   ierror )         &
       bind ( c, name='yac_cput_frac_' )

       use, intrinsic :: iso_c_binding, only : c_int, c_double

       integer ( kind=c_int ), value :: field_id
       integer ( kind=c_int ), value :: collection_size
       real    ( kind=c_double )     :: send_field(*)
       real    ( kind=c_double )     :: send_frac_mask(*)
       integer ( kind=c_int )        :: info
       integer ( kind=c_int )        :: ierror

     end subroutine yac_cput_frac__c

  end interface

  integer, intent (in)          :: field_id
  integer, intent (in)          :: nbr_hor_points
  integer, intent (in)          :: collection_size
  double precision, intent (in) :: send_field(nbr_hor_points, &
                                              collection_size)
  double precision, intent (in) :: send_frac_mask(nbr_hor_points, &
                                                  collection_size)
  integer, intent (out)         :: info
  integer, intent (out)         :: ierror

  call yac_fcheck_field_dimensions( &
    field_id, collection_size, 1, (/nbr_hor_points/) )
  YAC_FASSERT(is_contiguous(send_field), "ERROR(yac_fput_frac_single_pointset_dble): send_field is not contiguous")
  YAC_FASSERT(is_contiguous(send_frac_mask), "ERROR(yac_fput_frac_single_pointset_dble): send_frac_mask is not contiguous")

  call yac_cput_frac__c ( field_id,        &
                          collection_size, &
                          send_field,      &
                          send_frac_mask,  &
                          info,            &
                          ierror )

end subroutine yac_fput_frac_single_pointset_dble

! ---------------------------------- get -------------------------------

subroutine yac_fget_real ( field_id,        &
                           nbr_hor_points,  &
                           collection_size, &
                           recv_field,      &
                           info,            &
                           ierror )

  use mo_yac_finterface, dummy => yac_fget_real
  use mo_yac_real_to_dble_utils

  implicit none

  integer, intent (in)  :: field_id        !< [IN] field identifier
  integer, intent (in)  :: nbr_hor_points  !< [IN] number of horizontal points
  integer, intent (in)  :: collection_size !< [IN] number of vertical level or bundles
  real, intent (inout)  :: recv_field(nbr_hor_points, collection_size)
                                           !< [INOUT] returned field
  integer, intent (out) :: info            !< [OUT] returned info
  integer, intent (out) :: ierror          !< [OUT] returned error

  double precision :: recv_field_dble(nbr_hor_points, collection_size)

  call yac_fcheck_field_dimensions( &
    field_id, collection_size, 1, (/nbr_hor_points/) )

  call recv_field_to_dble(field_id,        &
                          nbr_hor_points,  &
                          collection_size, &
                          recv_field,      &
                          recv_field_dble)

  call yac_fget_dble ( field_id,        &
                       nbr_hor_points,  &
                       collection_size, &
                       recv_field_dble, &
                       info,            &
                       ierror )

  call recv_field_from_dble(field_id,        &
                            nbr_hor_points,  &
                            collection_size, &
                            recv_field_dble, &
                            recv_field)

end subroutine yac_fget_real

subroutine yac_fget_real_ptr ( field_id,        &
                               collection_size, &
                               recv_field,      &
                               info,            &
                               ierror )

  use mo_yac_finterface, dummy => yac_fget_real_ptr
  use mo_yac_real_to_dble_utils

  implicit none

  integer, intent (in)  :: field_id        !< [IN] field identifier
  integer, intent (in)  :: collection_size !< [IN] number of vertical level or bundles
  type(yac_real_ptr)    :: recv_field(collection_size) !< [INOUT] returned field
  integer, intent (out) :: info            !< [OUT] returned info
  integer, intent (out) :: ierror          !< [OUT] returned error handler

  type(yac_dble_ptr)    :: recv_field_dble(collection_size)

  call yac_fcheck_field_dimensions( &
    field_id, collection_size,  1, (/SIZE(recv_field(1)%p, 1)/) )

  call recv_field_to_dble_ptr(field_id,        &
                              collection_size, &
                              recv_field,      &
                              recv_field_dble)

  call yac_fget_dble_ptr ( field_id,        &
                           collection_size, &
                           recv_field_dble, &
                           info,            &
                           ierror )

  call recv_field_from_dble_ptr(field_id,        &
                                collection_size, &
                                recv_field_dble, &
                                recv_field)

end subroutine yac_fget_real_ptr

subroutine yac_fget_dble ( field_id,        &
                           nbr_hor_points,  &
                           collection_size, &
                           recv_field,      &
                           info,            &
                           ierror )

  use mo_yac_finterface, dummy => yac_fget_dble

  implicit none

  interface

     subroutine yac_cget__c ( field_id,        &
                              collection_size, &
                              recv_field,      &
                              info,            &
                              ierror ) bind ( c, name='yac_cget_' )

       use, intrinsic :: iso_c_binding, only : c_int, c_double

       integer ( kind=c_int ), value :: field_id            !< [IN] point identifier
       integer ( kind=c_int ), value :: collection_size     !< [IN] collection size
       real    ( kind=c_double )     :: recv_field(*)       !< [OUT] returned field
       integer ( kind=c_int )        :: info                !< [OUT] returned info
       integer ( kind=c_int )        :: ierror              !< [OUT] returned error handler

     end subroutine yac_cget__c

  end interface

  integer, intent (in)          :: field_id        !< [IN] field identifier
  integer, intent (in)          :: nbr_hor_points  !< [IN] number of horizontal points
  integer, intent (in)          :: collection_size !< [IN] number of vertical level or bundles
  double precision, intent (inout):: recv_field(nbr_hor_points, collection_size)
                                                  !< [INOUT] returned field
  integer, intent (out)         :: info           !< [OUT] returned info
  integer, intent (out)         :: ierror         !< [OUT] returned error

  call yac_fcheck_field_dimensions( &
    field_id, collection_size, 1, (/nbr_hor_points/) )
  YAC_FASSERT(is_contiguous(recv_field), "ERROR(yac_fget_dble): recv_field is not contiguous")

  call yac_cget__c ( field_id,        &
                     collection_size, &
                     recv_field,      &
                     info,            &
                     ierror )

end subroutine yac_fget_dble

subroutine yac_fget_dble_ptr ( field_id,        &
                               collection_size, &
                               recv_field,      &
                               info,            &
                               ierror )

  use mo_yac_finterface, dummy => yac_fget_dble_ptr
  use iso_c_binding, only: c_ptr, c_loc

  implicit none

  interface

     subroutine yac_cget_c ( field_id,        &
                             collection_size, &
                             recv_field,      &
                             info,            &
                             ierror ) bind ( c, name='yac_cget' )

       use, intrinsic :: iso_c_binding, only : c_int, c_ptr

       integer ( kind=c_int ), value :: field_id            !< [IN] point identifier
       integer ( kind=c_int ), value :: collection_size     !< [IN] collection size
       type(c_ptr)                   :: recv_field(*)       !< [OUT] returned field
       integer ( kind=c_int )        :: info                !< [OUT] returned info
       integer ( kind=c_int )        :: ierror              !< [OUT] returned error handler

     end subroutine yac_cget_c

  end interface

  integer, intent (in)  :: field_id        !< [IN] field identifier
  integer, intent (in)  :: collection_size !< [IN] number of vertical level or bundles
  type(yac_dble_ptr)      :: recv_field(collection_size) !< [OUT] returned field
  integer, intent (out) :: info            !< [OUT] returned info
  integer, intent (out) :: ierror          !< [OUT] returned error handler

  integer :: i
  type(c_ptr) :: recv_field_(collection_size)

  call yac_fcheck_field_dimensions( &
    field_id, collection_size, 1, (/SIZE(recv_field(1)%p)/) )

  do i = 1, collection_size
    recv_field_(i) = c_loc(recv_field(i)%p(1))
    YAC_FASSERT(is_contiguous(recv_field(i)%p), "ERROR(yac_fget_dble_ptr): recv_field is not contiguous")
  end do

  call yac_cget_c ( field_id,        &
                    collection_size, &
                    recv_field_,     &
                    info,            &
                    ierror )

end subroutine yac_fget_dble_ptr

! ---------------------------------- exchange --------------------------------

subroutine yac_fexchange_real ( send_field_id,       &
                                recv_field_id,       &
                                send_nbr_hor_points, &
                                send_nbr_pointsets,  &
                                recv_nbr_hor_points, &
                                collection_size,     &
                                send_field,          &
                                recv_field,          &
                                send_info,           &
                                recv_info,           &
                                ierror )

  use mo_yac_finterface, dummy => yac_fexchange_real
  use mo_yac_real_to_dble_utils

  implicit none

  integer, intent (in)  :: send_field_id       !< [IN] field identifier
  integer, intent (in)  :: recv_field_id       !< [IN] field identifier
  integer, intent (in)  :: send_nbr_hor_points !< [IN] number of horizontal send points
  integer, intent (in)  :: send_nbr_pointsets  !< [IN] number of send point sets
  integer, intent (in)  :: recv_nbr_hor_points !< [IN] number of horizontal recv points
  integer, intent (in)  :: collection_size     !< [IN] number of vertical level or bundles
  real, intent (in)     :: send_field(send_nbr_hor_points, &
                                      send_nbr_pointsets,  &
                                      collection_size)
                                               !< [IN] send field
  real, intent (inout)  :: recv_field(recv_nbr_hor_points, &
                                      collection_size)
                                               !< [INOUT] returned recv field
  integer, intent (out) :: send_info           !< [OUT] returned send info
  integer, intent (out) :: recv_info           !< [OUT] returned recv info
  integer, intent (out) :: ierror              !< [OUT] returned error

  double precision :: send_buffer(send_nbr_hor_points, &
                                  send_nbr_pointsets,  &
                                  collection_size)
  double precision :: recv_buffer(recv_nbr_hor_points, &
                                  collection_size)

  integer :: i

  call yac_fcheck_field_dimensions(                     &
    send_field_id, collection_size, send_nbr_pointsets, &
    (/(send_nbr_hor_points,i=1,send_nbr_pointsets)/) )
  call yac_fcheck_field_dimensions( &
    recv_field_id, collection_size, 1, (/recv_nbr_hor_points/) )

  call send_field_to_dble(send_field_id,       &
                          send_nbr_hor_points, &
                          send_nbr_pointsets,  &
                          collection_size,     &
                          send_field,          &
                          send_buffer)
  call recv_field_to_dble(recv_field_id,        &
                          recv_nbr_hor_points,  &
                          collection_size,      &
                          recv_field,           &
                          recv_buffer)

  call yac_fexchange_dble ( send_field_id,       &
                            recv_field_id,       &
                            send_nbr_hor_points, &
                            send_nbr_pointsets,  &
                            recv_nbr_hor_points, &
                            collection_size,     &
                            send_buffer,         &
                            recv_buffer,         &
                            send_info,           &
                            recv_info,           &
                            ierror )

  call recv_field_from_dble(recv_field_id,       &
                            recv_nbr_hor_points, &
                            collection_size,     &
                            recv_buffer,         &
                            recv_field)

end subroutine yac_fexchange_real

subroutine yac_fexchange_frac_real ( send_field_id,       &
                                     recv_field_id,       &
                                     send_nbr_hor_points, &
                                     send_nbr_pointsets,  &
                                     recv_nbr_hor_points, &
                                     collection_size,     &
                                     send_field,          &
                                     send_frac_mask,      &
                                     recv_field,          &
                                     send_info,           &
                                     recv_info,           &
                                     ierror )

  use mo_yac_finterface, dummy => yac_fexchange_frac_real
  use mo_yac_real_to_dble_utils

  implicit none

  integer, intent (in)  :: send_field_id       !< [IN] field identifier
  integer, intent (in)  :: recv_field_id       !< [IN] field identifier
  integer, intent (in)  :: send_nbr_hor_points !< [IN] number of horizontal send points
  integer, intent (in)  :: send_nbr_pointsets  !< [IN] number of send point sets
  integer, intent (in)  :: recv_nbr_hor_points !< [IN] number of horizontal recv points
  integer, intent (in)  :: collection_size     !< [IN] number of vertical level or bundles
  real, intent (in)     :: send_field(send_nbr_hor_points, &
                                      send_nbr_pointsets,  &
                                      collection_size)
                                               !< [IN] send field
  real, intent (in)     :: send_frac_mask(send_nbr_hor_points, &
                                          send_nbr_pointsets,  &
                                          collection_size)
                                               !< [IN] fractional mask
  real, intent (inout)  :: recv_field(recv_nbr_hor_points, &
                                      collection_size)
                                               !< [INOUT] returned recv field
  integer, intent (out) :: send_info           !< [OUT] returned send info
  integer, intent (out) :: recv_info           !< [OUT] returned recv info
  integer, intent (out) :: ierror              !< [OUT] returned error

  double precision :: send_buffer(send_nbr_hor_points, &
                                  send_nbr_pointsets,  &
                                  collection_size)
  double precision :: send_frac_mask_buffer(send_nbr_hor_points, &
                                            send_nbr_pointsets,  &
                                            collection_size)
  double precision :: recv_buffer(recv_nbr_hor_points, &
                                  collection_size)

  integer :: i

  call yac_fcheck_field_dimensions(                     &
    send_field_id, collection_size, send_nbr_pointsets, &
    (/(send_nbr_hor_points,i=1,send_nbr_pointsets)/) )
  call yac_fcheck_field_dimensions( &
    recv_field_id, collection_size, 1, (/recv_nbr_hor_points/) )

  call send_field_to_dble(send_field_id,       &
                          send_nbr_hor_points, &
                          send_nbr_pointsets,  &
                          collection_size,     &
                          send_field,          &
                          send_buffer,         &
                          send_frac_mask,      &
                          send_frac_mask_buffer)
  call recv_field_to_dble(recv_field_id,        &
                          recv_nbr_hor_points,  &
                          collection_size,      &
                          recv_field,           &
                          recv_buffer)

  call yac_fexchange_frac_dble ( send_field_id,         &
                                 recv_field_id,         &
                                 send_nbr_hor_points,   &
                                 send_nbr_pointsets,    &
                                 recv_nbr_hor_points,   &
                                 collection_size,       &
                                 send_buffer,           &
                                 send_frac_mask_buffer, &
                                 recv_buffer,           &
                                 send_info,             &
                                 recv_info,             &
                                 ierror )

  call recv_field_from_dble(recv_field_id,       &
                            recv_nbr_hor_points, &
                            collection_size,     &
                            recv_buffer,         &
                            recv_field)

end subroutine yac_fexchange_frac_real

subroutine yac_fexchange_real_ptr ( send_field_id,      &
                                    recv_field_id,      &
                                    send_nbr_pointsets, &
                                    collection_size,    &
                                    send_field,         &
                                    recv_field,         &
                                    send_info,          &
                                    recv_info,          &
                                    ierror )

  use mo_yac_finterface, dummy => yac_fexchange_real_ptr
  use mo_yac_real_to_dble_utils

  implicit none

  integer, intent (in)  :: send_field_id       !< [IN] send field identifier
  integer, intent (in)  :: recv_field_id       !< [IN] recv field identifier
  integer, intent (in)  :: send_nbr_pointsets  !< [IN] number of send point sets
  integer, intent (in)  :: collection_size     !< [IN] number of vertical level or bundles
  type(yac_real_ptr), intent (in) ::                      &
                           send_field(send_nbr_pointsets, &
                                      collection_size)
                                               !< [IN] send field
  type(yac_real_ptr)    :: recv_field(collection_size)
                                               !< [INOUT] returned recv field
  integer, intent (out) :: send_info           !< [OUT] returned send info
  integer, intent (out) :: recv_info           !< [OUT] returned recv info
  integer, intent (out) :: ierror              !< [OUT] returned error

  integer :: i, j
  type(yac_dble_ptr) :: send_field_dble(send_nbr_pointsets, collection_size)
  type(yac_dble_ptr) :: recv_field_dble(collection_size)

  call yac_fcheck_field_dimensions(                     &
    send_field_id, collection_size, send_nbr_pointsets, &
    (/(SIZE(send_field(i,1)%p),i=1,send_nbr_pointsets)/) )
  call yac_fcheck_field_dimensions( &
    recv_field_id, collection_size, 1, (/SIZE(recv_field(1)%p)/) )

  call send_field_to_dble_ptr(send_field_id,      &
                              send_nbr_pointsets, &
                              collection_size,    &
                              send_field,         &
                              send_field_dble)
  call recv_field_to_dble_ptr(recv_field_id,   &
                              collection_size, &
                              recv_field,      &
                              recv_field_dble)

  call yac_fexchange_dble_ptr ( send_field_id,      &
                                recv_field_id,      &
                                send_nbr_pointsets, &
                                collection_size,    &
                                send_field_dble,    &
                                recv_field_dble,    &
                                send_info,          &
                                recv_info,          &
                                ierror )

  call recv_field_from_dble_ptr(recv_field_id,   &
                                collection_size, &
                                recv_field_dble, &
                                recv_field)
  do i = 1, collection_size
    do j = 1, send_nbr_pointsets
      deallocate(send_field_dble(j, i)%p)
    end do
  end do

end subroutine yac_fexchange_real_ptr

subroutine yac_fexchange_frac_real_ptr ( send_field_id,      &
                                         recv_field_id,      &
                                         send_nbr_pointsets, &
                                         collection_size,    &
                                         send_field,         &
                                         send_frac_mask,     &
                                         recv_field,         &
                                         send_info,          &
                                         recv_info,          &
                                         ierror )

  use mo_yac_finterface, dummy => yac_fexchange_frac_real_ptr
  use mo_yac_real_to_dble_utils

  implicit none

  integer, intent (in)  :: send_field_id       !< [IN] send field identifier
  integer, intent (in)  :: recv_field_id       !< [IN] recv field identifier
  integer, intent (in)  :: send_nbr_pointsets  !< [IN] number of send point sets
  integer, intent (in)  :: collection_size     !< [IN] number of vertical level or bundles
  type(yac_real_ptr), intent (in) ::                      &
                           send_field(send_nbr_pointsets, &
                                      collection_size)
                                               !< [IN] send field
  type(yac_real_ptr), intent (in) ::                      &
                           send_frac_mask(send_nbr_pointsets, &
                                          collection_size)
                                               !< [IN] fractional mask
  type(yac_real_ptr)    :: recv_field(collection_size)
                                               !< [INOUT] returned recv field
  integer, intent (out) :: send_info           !< [OUT] returned send info
  integer, intent (out) :: recv_info           !< [OUT] returned recv info
  integer, intent (out) :: ierror              !< [OUT] returned error

  integer :: i, j
  type(yac_dble_ptr) :: send_field_dble(send_nbr_pointsets, collection_size)
  type(yac_dble_ptr) :: send_frac_mask_dble(send_nbr_pointsets, collection_size)
  type(yac_dble_ptr) :: recv_field_dble(collection_size)

  call yac_fcheck_field_dimensions(                     &
    send_field_id, collection_size, send_nbr_pointsets, &
    (/(SIZE(send_field(i,1)%p),i=1,send_nbr_pointsets)/) )
  call yac_fcheck_field_dimensions( &
    recv_field_id, collection_size, 1, (/SIZE(recv_field(1)%p)/) )

  call send_field_to_dble_ptr(send_field_id,      &
                              send_nbr_pointsets, &
                              collection_size,    &
                              send_field,         &
                              send_field_dble,    &
                              send_frac_mask,     &
                              send_frac_mask_dble)
  call recv_field_to_dble_ptr(recv_field_id,   &
                              collection_size, &
                              recv_field,      &
                              recv_field_dble)

  call yac_fexchange_frac_dble_ptr ( send_field_id,       &
                                     recv_field_id,       &
                                     send_nbr_pointsets,  &
                                     collection_size,     &
                                     send_field_dble,     &
                                     send_frac_mask_dble, &
                                     recv_field_dble,     &
                                     send_info,           &
                                     recv_info,           &
                                     ierror )

  call recv_field_from_dble_ptr(recv_field_id,   &
                                collection_size, &
                                recv_field_dble, &
                                recv_field)
  do i = 1, collection_size
    do j = 1, send_nbr_pointsets
      deallocate(send_field_dble(j, i)%p)
      deallocate(send_frac_mask_dble(j, i)%p)
    end do
  end do

end subroutine yac_fexchange_frac_real_ptr

subroutine yac_fexchange_single_pointset_real ( send_field_id,       &
                                                recv_field_id,       &
                                                send_nbr_hor_points, &
                                                recv_nbr_hor_points, &
                                                collection_size,     &
                                                send_field,          &
                                                recv_field,          &
                                                send_info,           &
                                                recv_info,           &
                                                ierror )

  use mo_yac_finterface, dummy => yac_fexchange_single_pointset_real
  use mo_yac_real_to_dble_utils

  implicit none

  integer, intent (in)  :: send_field_id        !< [IN] send field identifier
  integer, intent (in)  :: recv_field_id        !< [IN] recv field identifier
  integer, intent (in)  :: send_nbr_hor_points  !< [IN] number of horizontal send points
  integer, intent (in)  :: recv_nbr_hor_points  !< [IN] number of horizontal recv points
  integer, intent (in)  :: collection_size      !< [IN] number of vertical level or bundles
  real, intent (in)     :: send_field(send_nbr_hor_points, &
                                      collection_size)
                                                !< [IN] send field
  real, intent (inout)  :: recv_field(recv_nbr_hor_points, &
                                      collection_size)
                                                !< [INOUT] returned recv field
  integer, intent (out) :: send_info            !< [OUT] returned send info
  integer, intent (out) :: recv_info            !< [OUT] returned recv info
  integer, intent (out) :: ierror               !< [OUT] returned error

  double precision :: send_buffer(send_nbr_hor_points, collection_size)
  double precision :: recv_buffer(recv_nbr_hor_points, collection_size)

  call yac_fcheck_field_dimensions( &
    send_field_id, collection_size, 1, (/send_nbr_hor_points/) )
  call yac_fcheck_field_dimensions( &
    recv_field_id, collection_size, 1, (/recv_nbr_hor_points/) )

  call send_field_to_dble_single(send_field_id,        &
                                 send_nbr_hor_points,  &
                                 collection_size,      &
                                 send_field,           &
                                 send_buffer)

  call recv_field_to_dble(recv_field_id,       &
                          recv_nbr_hor_points, &
                          collection_size,     &
                          recv_field,          &
                          recv_buffer)

  call yac_fexchange_single_pointset_dble ( send_field_id,       &
                                            recv_field_id,       &
                                            send_nbr_hor_points, &
                                            recv_nbr_hor_points, &
                                            collection_size,     &
                                            send_buffer,         &
                                            recv_buffer,         &
                                            send_info,           &
                                            recv_info,           &
                                            ierror )

  call recv_field_from_dble(recv_field_id,       &
                            recv_nbr_hor_points, &
                            collection_size,     &
                            recv_buffer,         &
                            recv_field)

end subroutine yac_fexchange_single_pointset_real

subroutine yac_fexchange_frac_single_pointset_real ( send_field_id,       &
                                                     recv_field_id,       &
                                                     send_nbr_hor_points, &
                                                     recv_nbr_hor_points, &
                                                     collection_size,     &
                                                     send_field,          &
                                                     send_frac_mask,      &
                                                     recv_field,          &
                                                     send_info,           &
                                                     recv_info,           &
                                                     ierror )

  use mo_yac_finterface, dummy => yac_fexchange_frac_single_pointset_real
  use mo_yac_real_to_dble_utils

  implicit none

  integer, intent (in)  :: send_field_id        !< [IN] send field identifier
  integer, intent (in)  :: recv_field_id        !< [IN] recv field identifier
  integer, intent (in)  :: send_nbr_hor_points  !< [IN] number of horizontal send points
  integer, intent (in)  :: recv_nbr_hor_points  !< [IN] number of horizontal recv points
  integer, intent (in)  :: collection_size      !< [IN] number of vertical level or bundles
  real, intent (in)     :: send_field(send_nbr_hor_points, &
                                      collection_size)
                                                !< [IN] send field
  real, intent (in)     :: send_frac_mask(send_nbr_hor_points, &
                                          collection_size)
                                                !< [IN] fractional mask
  real, intent (inout)  :: recv_field(recv_nbr_hor_points, &
                                      collection_size)
                                                !< [INOUT] returned recv field
  integer, intent (out) :: send_info            !< [OUT] returned send info
  integer, intent (out) :: recv_info            !< [OUT] returned recv info
  integer, intent (out) :: ierror               !< [OUT] returned error

  double precision :: send_buffer(send_nbr_hor_points, collection_size)
  double precision :: send_frac_mask_buffer(send_nbr_hor_points, collection_size)
  double precision :: recv_buffer(recv_nbr_hor_points, collection_size)

  call yac_fcheck_field_dimensions( &
    send_field_id, collection_size, 1, (/send_nbr_hor_points/) )
  call yac_fcheck_field_dimensions( &
    recv_field_id, collection_size, 1, (/recv_nbr_hor_points/) )

  call send_field_to_dble_single(send_field_id,        &
                                 send_nbr_hor_points,  &
                                 collection_size,      &
                                 send_field,           &
                                 send_buffer,          &
                                 send_frac_mask,       &
                                 send_frac_mask_buffer)

  call recv_field_to_dble(recv_field_id,       &
                          recv_nbr_hor_points, &
                          collection_size,     &
                          recv_field,          &
                          recv_buffer)

  call yac_fexchange_frac_single_pointset_dble ( send_field_id,         &
                                                 recv_field_id,         &
                                                 send_nbr_hor_points,   &
                                                 recv_nbr_hor_points,   &
                                                 collection_size,       &
                                                 send_buffer,           &
                                                 send_frac_mask_buffer, &
                                                 recv_buffer,           &
                                                 send_info,             &
                                                 recv_info,             &
                                                 ierror )

  call recv_field_from_dble(recv_field_id,       &
                            recv_nbr_hor_points, &
                            collection_size,     &
                            recv_buffer,         &
                            recv_field)

end subroutine yac_fexchange_frac_single_pointset_real

!>
!! @param[in]    send_field_id       send field identifier
!! @param[in]    recv_field_id       recv field identifier
!! @param[in]    send_nbr_hor_points number of horizontal send points
!! @param[in]    send_nbr_pointsets  number of send point sets
!! @param[in]    recv_nbr_hor_points number of horizontal recv points
!! @param[in]    collection_size     number of vertical level or bundles
!! @param[in]    send_field          send field
!! @param[inout] recv_field          returned recv field
!! @param[out]   send_info           returned send info
!! @param[out]   recv_info           returned recv info
!! @param[out]   ierror              returned error
subroutine yac_fexchange_dble ( send_field_id,       &
                                recv_field_id,       &
                                send_nbr_hor_points, &
                                send_nbr_pointsets,  &
                                recv_nbr_hor_points, &
                                collection_size,     &
                                send_field,          &
                                recv_field,          &
                                send_info,           &
                                recv_info,           &
                                ierror )

  use mo_yac_finterface, dummy => yac_fexchange_dble

  implicit none

  interface

     subroutine yac_cexchange__c ( send_field_id,   &
                                   recv_field_id,   &
                                   collection_size, &
                                   send_field,      &
                                   recv_field,      &
                                   send_info,       &
                                   recv_info,       &
                                   ierror )         &
         bind ( c, name='yac_cexchange_' )

       use, intrinsic :: iso_c_binding, only : c_int, c_double

       integer ( kind=c_int ), value :: send_field_id
       integer ( kind=c_int ), value :: recv_field_id
       integer ( kind=c_int ), value :: collection_size
       real    ( kind=c_double )     :: send_field(*)
       real    ( kind=c_double )     :: recv_field(*)
       integer ( kind=c_int )        :: send_info
       integer ( kind=c_int )        :: recv_info
       integer ( kind=c_int )        :: ierror

     end subroutine yac_cexchange__c

  end interface

  integer, intent (in)  :: send_field_id
  integer, intent (in)  :: recv_field_id
  integer, intent (in)  :: send_nbr_hor_points
  integer, intent (in)  :: send_nbr_pointsets
  integer, intent (in)  :: recv_nbr_hor_points
  integer, intent (in)  :: collection_size
  double precision, intent (in) ::               &
                           send_field(           &
                            send_nbr_hor_points, &
                            send_nbr_pointsets,  &
                            collection_size)
  double precision, intent (inout)::             &
                           recv_field(           &
                            recv_nbr_hor_points, &
                            collection_size)
  integer, intent (out) :: send_info
  integer, intent (out) :: recv_info
  integer, intent (out) :: ierror

  integer :: i

  call yac_fcheck_field_dimensions( &
    send_field_id, collection_size, send_nbr_pointsets, &
    (/(send_nbr_hor_points,i=1,send_nbr_pointsets)/) )
  YAC_FASSERT(is_contiguous(send_field), "ERROR(yac_fexchange_dble): send_field is not contiguous")
  call yac_fcheck_field_dimensions( &
    recv_field_id, collection_size, 1, (/recv_nbr_hor_points/) )
  YAC_FASSERT(is_contiguous(recv_field), "ERROR(yac_fexchange_dble): recv_field is not contiguous")

  call yac_cexchange__c ( send_field_id,   &
                          recv_field_id,   &
                          collection_size, &
                          send_field,      &
                          recv_field,      &
                          send_info,       &
                          recv_info,       &
                          ierror )

end subroutine yac_fexchange_dble

!>
!! @param[in]    send_field_id       send field identifier
!! @param[in]    recv_field_id       recv field identifier
!! @param[in]    send_nbr_hor_points number of horizontal send points
!! @param[in]    send_nbr_pointsets  number of send point sets
!! @param[in]    recv_nbr_hor_points number of horizontal recv points
!! @param[in]    collection_size     number of vertical level or bundles
!! @param[in]    send_field          send field
!! @param[in]    send_frac_mask      fractional mask
!! @param[inout] recv_field          returned recv field
!! @param[out]   send_info           returned send info
!! @param[out]   recv_info           returned recv info
!! @param[out]   ierror              returned error
subroutine yac_fexchange_frac_dble ( send_field_id,       &
                                     recv_field_id,       &
                                     send_nbr_hor_points, &
                                     send_nbr_pointsets,  &
                                     recv_nbr_hor_points, &
                                     collection_size,     &
                                     send_field,          &
                                     send_frac_mask,      &
                                     recv_field,          &
                                     send_info,           &
                                     recv_info,           &
                                     ierror )

  use mo_yac_finterface, dummy => yac_fexchange_frac_dble

  implicit none

  interface

     subroutine yac_cexchange_frac__c ( send_field_id,   &
                                        recv_field_id,   &
                                        collection_size, &
                                        send_field,      &
                                        send_frac_mask,  &
                                        recv_field,      &
                                        send_info,       &
                                        recv_info,       &
                                        ierror )         &
         bind ( c, name='yac_cexchange_frac_' )

       use, intrinsic :: iso_c_binding, only : c_int, c_double

       integer ( kind=c_int ), value :: send_field_id
       integer ( kind=c_int ), value :: recv_field_id
       integer ( kind=c_int ), value :: collection_size
       real    ( kind=c_double )     :: send_field(*)
       real    ( kind=c_double )     :: send_frac_mask(*)
       real    ( kind=c_double )     :: recv_field(*)
       integer ( kind=c_int )        :: send_info
       integer ( kind=c_int )        :: recv_info
       integer ( kind=c_int )        :: ierror

     end subroutine yac_cexchange_frac__c

  end interface

  integer, intent (in)  :: send_field_id
  integer, intent (in)  :: recv_field_id
  integer, intent (in)  :: send_nbr_hor_points
  integer, intent (in)  :: send_nbr_pointsets
  integer, intent (in)  :: recv_nbr_hor_points
  integer, intent (in)  :: collection_size
  double precision, intent (in) ::               &
                           send_field(           &
                            send_nbr_hor_points, &
                            send_nbr_pointsets,  &
                            collection_size)
  double precision, intent (in) ::               &
                           send_frac_mask(           &
                            send_nbr_hor_points, &
                            send_nbr_pointsets,  &
                            collection_size)
  double precision, intent (inout)::             &
                           recv_field(           &
                            recv_nbr_hor_points, &
                            collection_size)
  integer, intent (out) :: send_info
  integer, intent (out) :: recv_info
  integer, intent (out) :: ierror

  integer :: i

  call yac_fcheck_field_dimensions( &
    send_field_id, collection_size, send_nbr_pointsets, &
    (/(send_nbr_hor_points,i=1,send_nbr_pointsets)/) )
  YAC_FASSERT(is_contiguous(send_field), "ERROR(yac_fexchange_frac_dble): send_field is not contiguous")
  YAC_FASSERT(is_contiguous(send_frac_mask), "ERROR(yac_fexchange_frac_dble): send_frac_mask is not contiguous")
  call yac_fcheck_field_dimensions( &
    recv_field_id, collection_size, 1, (/recv_nbr_hor_points/) )
  YAC_FASSERT(is_contiguous(recv_field), "ERROR(yac_fexchange_frac_dble): recv_field is not contiguous")

  call yac_cexchange_frac__c ( send_field_id,   &
                               recv_field_id,   &
                               collection_size, &
                               send_field,      &
                               send_frac_mask,  &
                               recv_field,      &
                               send_info,       &
                               recv_info,       &
                               ierror )

end subroutine yac_fexchange_frac_dble

!>
!! @param[in]    send_field_id       send field identifier
!! @param[in]    recv_field_id       recv field identifier
!! @param[in]    send_nbr_pointsets  number of send point sets
!! @param[in]    collection_size     number of vertical level or bundles
!! @param[in]    send_field          send field handle
!! @param[inout] recv_field          returned recv field handle
!! @param[out]   send_info           returned send info
!! @param[out]   recv_info           returned recv info
!! @param[out]   ierror              returned error
subroutine yac_fexchange_dble_ptr ( send_field_id,      &
                                    recv_field_id,      &
                                    send_nbr_pointsets, &
                                    collection_size,    &
                                    send_field,         &
                                    recv_field,         &
                                    send_info,          &
                                    recv_info,          &
                                    ierror )

  use mo_yac_finterface, dummy => yac_fexchange_dble_ptr
  use iso_c_binding, only: c_ptr, c_loc

  implicit none

  interface

     subroutine yac_cexchange_ptr__c ( send_field_id,   &
                                       recv_field_id,   &
                                       collection_size, &
                                       send_field,      &
                                       recv_field,      &
                                       send_info,       &
                                       recv_info,       &
                                       ierror )         &
         bind ( c, name='yac_cexchange_ptr_' )

       use, intrinsic :: iso_c_binding, only : c_int, c_ptr

       integer ( kind=c_int ), value :: send_field_id
       integer ( kind=c_int ), value :: recv_field_id
       integer ( kind=c_int ), value :: collection_size
       type(c_ptr)                   :: send_field(*)
       type(c_ptr)                   :: recv_field(*)
       integer ( kind=c_int )        :: send_info
       integer ( kind=c_int )        :: recv_info
       integer ( kind=c_int )        :: ierror

     end subroutine yac_cexchange_ptr__c

  end interface

  integer, intent (in)  :: send_field_id
  integer, intent (in)  :: recv_field_id
  integer, intent (in)  :: send_nbr_pointsets
  integer, intent (in)  :: collection_size
  type(yac_dble_ptr), intent (in) ::                      &
                           send_field(send_nbr_pointsets, &
                                      collection_size)
  type(yac_dble_ptr)      :: recv_field(collection_size)
  integer, intent (out) :: send_info
  integer, intent (out) :: recv_info
  integer, intent (out) :: ierror

  integer :: i, j
  type(c_ptr) :: send_field_(send_nbr_pointsets, collection_size)
  type(c_ptr) :: recv_field_(collection_size)

  call yac_fcheck_field_dimensions(                     &
    send_field_id, collection_size, send_nbr_pointsets, &
    (/(SIZE(send_field(i,1)%p),i=1,send_nbr_pointsets)/) )
  call yac_fcheck_field_dimensions( &
    recv_field_id, collection_size, 1, (/SIZE(recv_field(1)%p)/) )

  do i = 1, collection_size
    do j = 1, send_nbr_pointsets
      send_field_(j, i) = c_loc(send_field(j, i)%p(1))
      YAC_FASSERT(is_contiguous(send_field(j, i)%p), "ERROR(yac_fexchange_dble_ptr): send_field is not contiguous")
    end do
  end do
  do i = 1, collection_size
    recv_field_(i) = c_loc(recv_field(i)%p(1))
    YAC_FASSERT(is_contiguous(recv_field(i)%p), "ERROR(yac_fexchange_dble_ptr): recv_field is not contiguous")
  end do

  call yac_cexchange_ptr__c ( send_field_id,   &
                              recv_field_id,   &
                              collection_size, &
                              send_field_,     &
                              recv_field_,     &
                              send_info,       &
                              recv_info,       &
                              ierror )

end subroutine yac_fexchange_dble_ptr

!>
!! @param[in]    send_field_id       send field identifier
!! @param[in]    recv_field_id       recv field identifier
!! @param[in]    send_nbr_pointsets  number of send point sets
!! @param[in]    collection_size     number of vertical level or bundles
!! @param[in]    send_field          send field handle
!! @param[in]    send_frac_mask      fractional mask handle
!! @param[inout] recv_field          returned recv field handle
!! @param[out]   send_info           returned send info
!! @param[out]   recv_info           returned recv info
!! @param[out]   ierror              returned error
subroutine yac_fexchange_frac_dble_ptr ( send_field_id,      &
                                         recv_field_id,      &
                                         send_nbr_pointsets, &
                                         collection_size,    &
                                         send_field,         &
                                         send_frac_mask,     &
                                         recv_field,         &
                                         send_info,          &
                                         recv_info,          &
                                         ierror )

  use mo_yac_finterface, dummy => yac_fexchange_frac_dble_ptr
  use iso_c_binding, only: c_ptr, c_loc

  implicit none

  interface

     subroutine yac_cexchange_frac_ptr__c ( send_field_id,   &
                                            recv_field_id,   &
                                            collection_size, &
                                            send_field,      &
                                            send_frac_mask,  &
                                            recv_field,      &
                                            send_info,       &
                                            recv_info,       &
                                            ierror )         &
         bind ( c, name='yac_cexchange_frac_ptr_' )

       use, intrinsic :: iso_c_binding, only : c_int, c_ptr

       integer ( kind=c_int ), value :: send_field_id
       integer ( kind=c_int ), value :: recv_field_id
       integer ( kind=c_int ), value :: collection_size
       type(c_ptr)                   :: send_field(*)
       type(c_ptr)                   :: send_frac_mask(*)
       type(c_ptr)                   :: recv_field(*)
       integer ( kind=c_int )        :: send_info
       integer ( kind=c_int )        :: recv_info
       integer ( kind=c_int )        :: ierror

     end subroutine yac_cexchange_frac_ptr__c

  end interface

  integer, intent (in)  :: send_field_id
  integer, intent (in)  :: recv_field_id
  integer, intent (in)  :: send_nbr_pointsets
  integer, intent (in)  :: collection_size
  type(yac_dble_ptr), intent (in) ::                      &
                           send_field(send_nbr_pointsets, &
                                      collection_size)
  type(yac_dble_ptr), intent (in) ::                      &
                           send_frac_mask(send_nbr_pointsets, &
                                          collection_size)
  type(yac_dble_ptr)      :: recv_field(collection_size)
  integer, intent (out) :: send_info
  integer, intent (out) :: recv_info
  integer, intent (out) :: ierror

  integer :: i, j
  type(c_ptr) :: send_field_(send_nbr_pointsets, collection_size)
  type(c_ptr) :: send_frac_mask_(send_nbr_pointsets, collection_size)
  type(c_ptr) :: recv_field_(collection_size)

  call yac_fcheck_field_dimensions(                     &
    send_field_id, collection_size, send_nbr_pointsets, &
    (/(SIZE(send_field(i,1)%p),i=1,send_nbr_pointsets)/) )
  call yac_fcheck_field_dimensions(                     &
    send_field_id, collection_size, send_nbr_pointsets, &
    (/(SIZE(send_frac_mask(i,1)%p),i=1,send_nbr_pointsets)/) )
  call yac_fcheck_field_dimensions( &
    recv_field_id, collection_size, 1, (/SIZE(recv_field(1)%p)/) )

  do i = 1, collection_size
    do j = 1, send_nbr_pointsets
      send_field_(j, i) = c_loc(send_field(j, i)%p(1))
      YAC_FASSERT(is_contiguous(send_field(j, i)%p), "ERROR(yac_fexchange_frac_dble_ptr): send_field is not contiguous")
      send_frac_mask_(j, i) = c_loc(send_frac_mask(j, i)%p(1))
      YAC_FASSERT(is_contiguous(send_frac_mask(j, i)%p), "ERROR(yac_fexchange_frac_dble_ptr): send_frac_mask is not contiguous")
    end do
  end do
  do i = 1, collection_size
    recv_field_(i) = c_loc(recv_field(i)%p(1))
    YAC_FASSERT(is_contiguous(recv_field(i)%p), "ERROR(yac_fexchange_frac_dble_ptr): recv_field is not contiguous")
  end do

  call yac_cexchange_frac_ptr__c ( send_field_id,   &
                                   recv_field_id,   &
                                   collection_size, &
                                   send_field_,     &
                                   send_frac_mask_, &
                                   recv_field_,     &
                                   send_info,       &
                                   recv_info,       &
                                   ierror )

end subroutine yac_fexchange_frac_dble_ptr

!>
!! @param[in]    send_field_id       send field identifier
!! @param[in]    recv_field_id       recv field identifier
!! @param[in]    send_nbr_hor_points number of horizontal send points
!! @param[in]    recv_nbr_hor_points number of horizontal recv points
!! @param[in]    collection_size     number of vertical level or bundles
!! @param[in]    send_field          send field
!! @param[inout] recv_field          returned recv field
!! @param[out]   send_info           returned send info
!! @param[out]   recv_info           returned recv info
!! @param[out]   ierror              returned error
subroutine yac_fexchange_single_pointset_dble ( send_field_id,       &
                                                recv_field_id,       &
                                                send_nbr_hor_points, &
                                                recv_nbr_hor_points, &
                                                collection_size,     &
                                                send_field,          &
                                                recv_field,          &
                                                send_info,           &
                                                recv_info,           &
                                                ierror )

  use mo_yac_finterface, dummy => yac_fexchange_single_pointset_dble

  implicit none

  interface

     subroutine yac_cexchange__c ( send_field_id,   &
                                   recv_field_id,   &
                                   collection_size, &
                                   send_field,      &
                                   recv_field,      &
                                   send_info,       &
                                   recv_info,       &
                                   ierror )         &
         bind ( c, name='yac_cexchange_' )

       use, intrinsic :: iso_c_binding, only : c_int, c_double

       integer ( kind=c_int ), value :: send_field_id
       integer ( kind=c_int ), value :: recv_field_id
       integer ( kind=c_int ), value :: collection_size
       real    ( kind=c_double )     :: send_field(*)
       real    ( kind=c_double )     :: recv_field(*)
       integer ( kind=c_int )        :: send_info
       integer ( kind=c_int )        :: recv_info
       integer ( kind=c_int )        :: ierror

     end subroutine yac_cexchange__c

  end interface

  integer, intent (in)  :: send_field_id
  integer, intent (in)  :: recv_field_id
  integer, intent (in)  :: send_nbr_hor_points
  integer, intent (in)  :: recv_nbr_hor_points
  integer, intent (in)  :: collection_size
  double precision, intent (in) ::                         &
                           send_field(send_nbr_hor_points, &
                                      collection_size)
  double precision, intent (inout)::                       &
                           recv_field(recv_nbr_hor_points, &
                                      collection_size)
  integer, intent (out) :: send_info
  integer, intent (out) :: recv_info
  integer, intent (out) :: ierror

  call yac_fcheck_field_dimensions( &
    send_field_id, collection_size, 1, (/send_nbr_hor_points/) )
  YAC_FASSERT(is_contiguous(send_field), "ERROR(yac_fexchange_single_pointset_dble): send_field is not contiguous")
  call yac_fcheck_field_dimensions( &
    recv_field_id, collection_size, 1, (/recv_nbr_hor_points/) )
  YAC_FASSERT(is_contiguous(recv_field), "ERROR(yac_fexchange_single_pointset_dble): recv_field is not contiguous")

  call yac_cexchange__c ( send_field_id,   &
                          recv_field_id,   &
                          collection_size, &
                          send_field,      &
                          recv_field,      &
                          send_info,       &
                          recv_info,       &
                          ierror )

end subroutine yac_fexchange_single_pointset_dble

!>
!! @param[in]    send_field_id       send field identifier
!! @param[in]    recv_field_id       recv field identifier
!! @param[in]    send_nbr_hor_points number of horizontal send points
!! @param[in]    recv_nbr_hor_points number of horizontal recv points
!! @param[in]    collection_size     number of vertical level or bundles
!! @param[in]    send_field          send field
!! @param[in]    send_frac_mask      fractional mask
!! @param[inout] recv_field          returned recv field
!! @param[out]   send_info           returned send info
!! @param[out]   recv_info           returned recv info
!! @param[out]   ierror              returned error
subroutine yac_fexchange_frac_single_pointset_dble ( send_field_id,       &
                                                     recv_field_id,       &
                                                     send_nbr_hor_points, &
                                                     recv_nbr_hor_points, &
                                                     collection_size,     &
                                                     send_field,          &
                                                     send_frac_mask,      &
                                                     recv_field,          &
                                                     send_info,           &
                                                     recv_info,           &
                                                     ierror )

  use mo_yac_finterface, dummy => yac_fexchange_frac_single_pointset_dble

  implicit none

  interface

     subroutine yac_cexchange_frac__c ( send_field_id,   &
                                        recv_field_id,   &
                                        collection_size, &
                                        send_field,      &
                                        send_frac_mask,  &
                                        recv_field,      &
                                        send_info,       &
                                        recv_info,       &
                                        ierror )         &
         bind ( c, name='yac_cexchange_frac_' )

       use, intrinsic :: iso_c_binding, only : c_int, c_double

       integer ( kind=c_int ), value :: send_field_id
       integer ( kind=c_int ), value :: recv_field_id
       integer ( kind=c_int ), value :: collection_size
       real    ( kind=c_double )     :: send_field(*)
       real    ( kind=c_double )     :: send_frac_mask(*)
       real    ( kind=c_double )     :: recv_field(*)
       integer ( kind=c_int )        :: send_info
       integer ( kind=c_int )        :: recv_info
       integer ( kind=c_int )        :: ierror

     end subroutine yac_cexchange_frac__c

  end interface

  integer, intent (in)  :: send_field_id
  integer, intent (in)  :: recv_field_id
  integer, intent (in)  :: send_nbr_hor_points
  integer, intent (in)  :: recv_nbr_hor_points
  integer, intent (in)  :: collection_size
  double precision, intent (in) ::                         &
                           send_field(send_nbr_hor_points, &
                                      collection_size)
  double precision, intent (in) ::                         &
                           send_frac_mask(send_nbr_hor_points, &
                                          collection_size)
  double precision, intent (inout)::                       &
                           recv_field(recv_nbr_hor_points, &
                                      collection_size)
  integer, intent (out) :: send_info
  integer, intent (out) :: recv_info
  integer, intent (out) :: ierror

  call yac_fcheck_field_dimensions( &
    send_field_id, collection_size, 1, (/send_nbr_hor_points/) )
  YAC_FASSERT(is_contiguous(send_field), "ERROR(yac_fexchange_frac_single_pointset_dble): send_field is not contiguous")
  YAC_FASSERT(is_contiguous(send_frac_mask), "ERROR(yac_fexchange_frac_single_pointset_dble): send_frac_mask is not contiguous")
  call yac_fcheck_field_dimensions( &
    recv_field_id, collection_size, 1, (/recv_nbr_hor_points/) )
  YAC_FASSERT(is_contiguous(recv_field), "ERROR(yac_fexchange_frac_single_pointset_dble): recv_field is not contiguous")

  call yac_cexchange_frac__c ( send_field_id,   &
                               recv_field_id,   &
                               collection_size, &
                               send_field,      &
                               send_frac_mask,  &
                               recv_field,      &
                               send_info,       &
                               recv_info,       &
                               ierror )

end subroutine yac_fexchange_frac_single_pointset_dble

! ----------------------------------------------------------------------

subroutine yac_ftest_i ( field_id, flag )

  use mo_yac_finterface, dummy => yac_ftest_i

  implicit none

  interface

     subroutine yac_ctest_c ( field_id, flag ) &
         bind ( c, name='yac_ctest' )

       use, intrinsic :: iso_c_binding, only : c_int

       integer ( kind=c_int ), value :: field_id
       integer ( kind=c_int )        :: flag

     end subroutine yac_ctest_c

  end interface

  integer, intent (in)  :: field_id
  integer, intent (out) :: flag

  call yac_ctest_c ( field_id, flag )

end subroutine yac_ftest_i

subroutine yac_ftest_l ( field_id, flag )

  use mo_yac_finterface, dummy => yac_ftest_l

  implicit none

  interface

     subroutine yac_ctest_c ( field_id, flag ) &
         bind ( c, name='yac_ctest' )

       use, intrinsic :: iso_c_binding, only : c_int

       integer ( kind=c_int ), value :: field_id
       integer ( kind=c_int )        :: flag

     end subroutine yac_ctest_c

  end interface

  integer, intent (in)  :: field_id
  logical, intent (out) :: flag

  integer :: iflag

  call yac_ctest_c ( field_id, iflag )

  flag = iflag /= 0

end subroutine yac_ftest_l

! ----------------------------------------------------------------------

subroutine yac_fget_comp_comm ( comp_id, comp_comm )

   use mo_yac_finterface, dummy => yac_fget_comp_comm

   implicit none

   interface

     subroutine yac_get_comp_comm_c ( comp_id, comp_comm ) &
       bind ( c, name='yac_get_comp_comm_f2c' )

       use, intrinsic :: iso_c_binding, only : c_int
       use mo_yac_finterface, only : YAC_MPI_FINT_KIND

       integer ( kind=c_int ), value      :: comp_id
       integer ( kind=YAC_MPI_FINT_KIND ) :: comp_comm

     end subroutine yac_get_comp_comm_c

   end interface

   integer, intent (in)  :: comp_id   !< [IN] component identifier
   integer, intent (out) :: comp_comm !< [OUT] component MPI communicator

   call yac_get_comp_comm_c ( comp_id, comp_comm )

end subroutine yac_fget_comp_comm

! ----------------------------------------------------------------------

subroutine yac_fget_comps_comm ( comp_names, num_comps, comps_comm )

   use, intrinsic :: iso_c_binding, only : c_null_char, c_ptr, c_loc, c_char
   use mo_yac_finterface, dummy => yac_fget_comps_comm
   use mo_yac_iso_c_helpers

   implicit none

   interface

     subroutine yac_cget_comps_comm_c ( comp_names, &
                                        num_comps,  &
                                        comps_comm) &
           bind ( c, name='yac_cget_comps_comm_f2c' )

       use, intrinsic :: iso_c_binding, only : c_int, c_ptr
       use mo_yac_finterface, only : YAC_MPI_FINT_KIND

       type ( c_ptr )                     :: comp_names(*)
       integer ( kind=c_int ), value      :: num_comps
       integer ( kind=YAC_MPI_FINT_KIND ) :: comps_comm

     end subroutine yac_cget_comps_comm_c

   end interface

   integer, intent(in)   :: num_comps             !< [IN]  number of components
   character(kind=c_char, len=*), intent(in) :: &
                            comp_names(num_comps) !< [IN]  component names
   integer, intent (out) :: comps_comm            !< [OUT] components MPI communicator

   integer :: i, j
   character(kind=c_char), target :: comp_names_cpy(YAC_MAX_CHARLEN+1, num_comps)
   type(c_ptr) :: comp_name_ptrs(num_comps)

   comp_names_cpy = c_null_char

   do i = 1, num_comps
     call yac_check_strlength(comp_names(i))
     do j = 1, len_trim(comp_names(i))
       comp_names_cpy(j,i) = comp_names(i)(j:j)
     end do
     comp_name_ptrs(i) = c_loc(comp_names_cpy(1,i))
   end do

   call yac_cget_comps_comm_c ( comp_name_ptrs, &
                                num_comps,      &
                                comps_comm )

end subroutine yac_fget_comps_comm

subroutine yac_fget_comps_comm_instance ( yac_instance_id, &
                                          comp_names,      &
                                          num_comps,       &
                                          comps_comm )

   use, intrinsic :: iso_c_binding, only : c_null_char, c_ptr, c_loc, c_char
   use mo_yac_finterface, dummy => yac_fget_comps_comm_instance
   use mo_yac_iso_c_helpers

   implicit none

   interface

     subroutine yac_cget_comps_comm_instance_c ( yac_instance_id, &
                                                 comp_names,      &
                                                 num_comps,       &
                                                 comps_comm)      &
           bind ( c, name='yac_cget_comps_comm_instance_f2c' )

       use, intrinsic :: iso_c_binding, only : c_int, c_ptr
       use mo_yac_finterface, only : YAC_MPI_FINT_KIND

       integer   ( kind=c_int ), value    :: yac_instance_id
       type ( c_ptr )                     :: comp_names(*)
       integer ( kind=c_int ), value      :: num_comps
       integer ( kind=YAC_MPI_FINT_KIND ) :: comps_comm

     end subroutine yac_cget_comps_comm_instance_c

   end interface

   integer, intent(in)   :: yac_instance_id       !< [IN]  YAC instance identifier
   integer, intent(in)   :: num_comps             !< [IN]  number of components
   character(kind=c_char, len=*), intent(in) :: &
                            comp_names(num_comps) !< [IN]  component names
   integer, intent (out) :: comps_comm            !< [OUT] components MPI communicator

   integer :: i, j
   character(kind=c_char), target :: comp_names_cpy(YAC_MAX_CHARLEN+1, num_comps)
   type(c_ptr) :: comp_name_ptrs(num_comps)

   comp_names_cpy = c_null_char

   do i = 1, num_comps
     call yac_check_strlength(comp_names(i))
     do j = 1, len_trim(comp_names(i))
       comp_names_cpy(j,i) = comp_names(i)(j:j)
     end do
     comp_name_ptrs(i) = c_loc(comp_names_cpy(1,i))
   end do

   call yac_cget_comps_comm_instance_c ( yac_instance_id, &
                                         comp_name_ptrs,  &
                                         num_comps,       &
                                         comps_comm )

end subroutine yac_fget_comps_comm_instance

! ------------------- search/end of definition -------------------------

subroutine yac_fsync_def ( )

   use mo_yac_finterface, dummy => yac_fsync_def

   implicit none

   interface

     subroutine yac_csync_def_c ( ) bind ( c, name='yac_csync_def' )

     end subroutine yac_csync_def_c

   end interface

   call yac_csync_def_c ( )

end subroutine yac_fsync_def

subroutine yac_fsync_def_instance ( yac_instance_id )

   use mo_yac_finterface, dummy => yac_fsync_def_instance

   implicit none

   interface

     subroutine yac_csync_def_instance_c ( yac_instance_id ) &
       bind ( c, name='yac_csync_def_instance' )

       use, intrinsic :: iso_c_binding, only : c_int

       integer ( kind=c_int ), value :: yac_instance_id

     end subroutine yac_csync_def_instance_c

   end interface

   integer, intent(in)  :: yac_instance_id !< [IN]  YAC instance identifier

   call yac_csync_def_instance_c ( yac_instance_id )

end subroutine yac_fsync_def_instance


subroutine yac_fenddef ( )

   use mo_yac_finterface, dummy => yac_fenddef

   implicit none

   interface

     subroutine yac_cenddef_c ( ) bind ( c, name='yac_cenddef' )

     end subroutine yac_cenddef_c

   end interface

   call yac_cenddef_c ( )

end subroutine yac_fenddef

subroutine yac_fenddef_instance ( yac_instance_id )

   use mo_yac_finterface, dummy => yac_fenddef_instance

   implicit none

   interface

     subroutine yac_cenddef_instance_c ( yac_instance_id ) &
       bind ( c, name='yac_cenddef_instance' )

       use, intrinsic :: iso_c_binding, only : c_int

       integer ( kind=c_int ), value :: yac_instance_id

     end subroutine yac_cenddef_instance_c

   end interface

   integer, intent(in)  :: yac_instance_id !< [IN]  YAC instance identifier

   call yac_cenddef_instance_c ( yac_instance_id )

end subroutine yac_fenddef_instance

subroutine yac_fenddef_and_emit_config(emit_flags, config)

  use, intrinsic :: iso_c_binding, only : c_ptr
  use mo_yac_finterface, dummy => yac_fenddef_and_emit_config
  use mo_yac_iso_c_helpers

  implicit none

  interface
     subroutine yac_fenddef_and_emit_config_c ( &
       emit_flags, config)                    &
       bind ( c, name='yac_cenddef_and_emit_config' )

     use, intrinsic :: iso_c_binding, only : c_int, c_ptr
     integer ( kind=c_int ), value :: emit_flags
     type(c_ptr) :: config

     end subroutine yac_fenddef_and_emit_config_c

     subroutine free_c ( ptr ) BIND ( c, NAME='free' )

       use, intrinsic :: iso_c_binding, only : c_ptr

       type ( c_ptr ), intent(in), value :: ptr

     end subroutine free_c
  end interface

  integer, intent (in)           :: emit_flags !< [IN] flags for emitting the config
  character (len=:), ALLOCATABLE :: config     !< [IN,OUT] configuration string

  type (c_ptr)                   :: c_string_ptr

  call yac_fenddef_and_emit_config_c(emit_flags, c_string_ptr)
  config = yac_internal_cptr2char(c_string_ptr)
  call free_c(c_string_ptr)

end subroutine yac_fenddef_and_emit_config

subroutine yac_fenddef_and_emit_config_instance( &
  yac_instance_id, emit_flags, config)

  use, intrinsic :: iso_c_binding, only : c_ptr
  use mo_yac_finterface, dummy => yac_fenddef_and_emit_config_instance
  use mo_yac_iso_c_helpers

  implicit none

  interface
     subroutine yac_fenddef_and_emit_config_instance_c ( &
       yac_instance_id, emit_flags, config)            &
       bind ( c, name='yac_cenddef_and_emit_config_instance' )

     use, intrinsic :: iso_c_binding, only : c_int, c_ptr
     integer ( kind=c_int ), value :: yac_instance_id
     integer ( kind=c_int ), value :: emit_flags
     type(c_ptr) :: config

     end subroutine yac_fenddef_and_emit_config_instance_c

     subroutine free_c ( ptr ) BIND ( c, NAME='free' )

       use, intrinsic :: iso_c_binding, only : c_ptr

       type ( c_ptr ), intent(in), value :: ptr

     end subroutine free_c
  end interface

  integer, intent (in)           :: yac_instance_id !< [IN] YAC instance identifier
  integer, intent (in)           :: emit_flags      !< [IN] flags for emitting the config
  character (len=:), ALLOCATABLE :: config          !< [IN,OUT] configuration string

  type (c_ptr)                   :: c_string_ptr

  call yac_fenddef_and_emit_config_instance_c( &
    yac_instance_id, emit_flags, c_string_ptr)
  config = yac_internal_cptr2char(c_string_ptr)
  call free_c(c_string_ptr)

end subroutine yac_fenddef_and_emit_config_instance

! ------------------------ query routines -----------------------------

 function yac_fget_comp_names ( ) result( comp_names )

   use mo_yac_finterface, dummy => yac_fget_comp_names
   use mo_yac_iso_c_helpers

   use, intrinsic :: iso_c_binding, ONLY: c_ptr
   implicit none

   interface
      function yac_cget_nbr_comps_c() result( nbr_comps ) &
           bind(c, name='yac_cget_nbr_comps')
        use, intrinsic :: iso_c_binding, only: c_int
        integer(kind=c_int) :: nbr_comps

      end function yac_cget_nbr_comps_c

      subroutine yac_cget_comp_names_c( nbr_comps, comp_names ) &
           bind(c, name='yac_cget_comp_names')
        use, intrinsic :: iso_c_binding, ONLY: c_int, c_ptr
        integer(kind=c_int), intent(in), value :: nbr_comps
        TYPE(c_ptr), intent(out) :: comp_names(nbr_comps)
      end subroutine yac_cget_comp_names_c
   end interface

   type(yac_string), allocatable :: comp_names(:)
   integer :: nbr_comps
   INTEGER :: i
   TYPE(c_ptr), allocatable :: comp_ptr(:)

   nbr_comps = yac_cget_nbr_comps_c()
   allocate(comp_ptr(nbr_comps))
   allocate(comp_names(nbr_comps))
   CALL yac_cget_comp_names_c(nbr_comps, comp_ptr)
   DO i=1,nbr_comps
      comp_names(i)%string = yac_internal_cptr2char(comp_ptr(i))
   END DO
 end function yac_fget_comp_names

 function yac_fget_comp_names_instance ( yac_instance_id ) result ( comp_names )

   use mo_yac_finterface, dummy => yac_fget_comp_names_instance
   use mo_yac_iso_c_helpers
   use, intrinsic :: iso_c_binding, ONLY: c_ptr
   implicit none

   interface
      function yac_cget_nbr_comps_instance_c( yac_instance_id ) &
          result( nbr_comps )                                   &
           bind(c, name='yac_cget_nbr_comps_instance')
        use, intrinsic :: iso_c_binding, only: c_int
        integer(kind=c_int), value, intent(in) :: yac_instance_id
        integer(kind=c_int) :: nbr_comps
      end function yac_cget_nbr_comps_instance_c

      subroutine yac_cget_comp_names_instance_c( yac_instance_id, &
                                                 nbr_comps,       &
                                                 comp_names )     &
           bind(c, name='yac_cget_comp_names_instance')
        use, intrinsic :: iso_c_binding, ONLY: c_int, c_ptr
        integer(kind=c_int), intent(in), value :: yac_instance_id
        integer(kind=c_int), intent(in), value :: nbr_comps
        TYPE(c_ptr), intent(out) :: comp_names(nbr_comps)
      end subroutine yac_cget_comp_names_instance_c
   end interface

   integer, intent(in) :: yac_instance_id
   type(yac_string), allocatable :: comp_names(:)
   integer :: nbr_comps
   INTEGER :: i
   TYPE(c_ptr), allocatable :: comp_ptr(:)

   nbr_comps = yac_cget_nbr_comps_instance_c(yac_instance_id)
   allocate(comp_names(nbr_comps))
   allocate(comp_ptr(nbr_comps))
   CALL yac_cget_comp_names_instance_c(yac_instance_id, &
                                       nbr_comps,       &
                                       comp_ptr)
   DO i=1,nbr_comps
      comp_names(i)%string = yac_internal_cptr2char(comp_ptr(i))
   END DO
 end function  yac_fget_comp_names_instance

! ---------------------------------------------------------------------

 function yac_fget_grid_names ( ) result ( grid_names )

   use mo_yac_finterface, dummy => yac_fget_grid_names
   use mo_yac_iso_c_helpers
   use, intrinsic :: iso_c_binding, ONLY: c_ptr

   implicit none

   interface
      function yac_cget_nbr_grids_c() result( nbr_grids ) &
           bind(c, name='yac_cget_nbr_grids')
        use, intrinsic :: iso_c_binding, only: c_int
        integer(kind=c_int) :: nbr_grids
      end function yac_cget_nbr_grids_c

      subroutine yac_cget_grid_names_c( nbr_grids, grid_names ) &
           bind(c, name='yac_cget_grid_names')
        use, intrinsic :: iso_c_binding, ONLY: c_int, c_ptr
        integer(kind=c_int), intent(in), value :: nbr_grids
        TYPE(c_ptr), intent(out) :: grid_names(nbr_grids)
      end subroutine yac_cget_grid_names_c
   end interface

   type(yac_string), allocatable :: grid_names(:)
   integer :: nbr_grids
   INTEGER :: i
   TYPE(c_ptr), allocatable :: grid_ptr(:)

   nbr_grids = yac_cget_nbr_grids_c()
   allocate(grid_ptr(nbr_grids))
   CALL yac_cget_grid_names_c(nbr_grids, grid_ptr)
   allocate(grid_names(nbr_grids))
   DO i=1,nbr_grids
      grid_names(i)%string = yac_internal_cptr2char(grid_ptr(i))
   END DO
 end function yac_fget_grid_names

 function yac_fget_grid_names_instance ( yac_instance_id ) result ( grid_names )

   use mo_yac_finterface, dummy => yac_fget_grid_names_instance
   use mo_yac_iso_c_helpers
   use, intrinsic :: iso_c_binding, ONLY: c_ptr

   implicit none

   interface
      function yac_cget_nbr_grids_instance_c( yac_instance_id ) &
           result( nbr_grids )                                  &
           bind(c, name='yac_cget_nbr_grids_instance')
        use, intrinsic :: iso_c_binding, only: c_int
        integer(kind=c_int), value, intent(in) :: yac_instance_id
        integer(kind=c_int) :: nbr_grids

      end function yac_cget_nbr_grids_instance_c

      subroutine yac_cget_grid_names_instance_c( yac_instance_id, &
                                                 nbr_grids,       &
                                                 grid_names )     &
           bind(c, name='yac_cget_grid_names_instance')
        use, intrinsic :: iso_c_binding, ONLY: c_int, c_ptr
        integer(kind=c_int), intent(in), value :: yac_instance_id
        integer(kind=c_int), intent(in), value :: nbr_grids
        TYPE(c_ptr), intent(out) :: grid_names(nbr_grids)
      end subroutine yac_cget_grid_names_instance_c
   end interface

   integer, intent(in) :: yac_instance_id
   type(yac_string), allocatable :: grid_names(:)
   integer :: nbr_grids
   INTEGER :: i
   TYPE(c_ptr), allocatable :: grid_ptr(:)

   nbr_grids = yac_cget_nbr_grids_instance_c(yac_instance_id)
   allocate(grid_ptr(nbr_grids))
   CALL yac_cget_grid_names_instance_c(yac_instance_id, &
                                       nbr_grids,       &
                                       grid_ptr)
   allocate(grid_names(nbr_grids))
   DO i=1,nbr_grids
      grid_names(i)%string = yac_internal_cptr2char(grid_ptr(i))
   END DO
 end function yac_fget_grid_names_instance

 ! ---------------------------------------------------------------------

 function yac_fget_comp_grid_names ( comp_name ) result ( grid_names )

   use mo_yac_finterface, dummy => yac_fget_comp_grid_names
   use mo_yac_iso_c_helpers
   use, intrinsic :: iso_c_binding, ONLY: c_ptr, c_null_char

   implicit none

   interface
      function yac_cget_comp_nbr_grids_c( comp_name ) result( nbr_grids ) &
           bind(c, name='yac_cget_comp_nbr_grids')
        use, intrinsic :: iso_c_binding, only: c_int, c_char
        character(kind=c_char), dimension(*), intent(in) :: comp_name
        integer(kind=c_int) :: nbr_grids

      end function yac_cget_comp_nbr_grids_c

      subroutine yac_cget_comp_grid_names_c( comp_name,   &
                                             nbr_grids,   &
                                             grid_names ) &
           bind(c, name='yac_cget_comp_grid_names')
        use, intrinsic :: iso_c_binding, ONLY: c_int, c_ptr, c_char
        character(kind=c_char), dimension(*), intent(in) :: comp_name
        integer(kind=c_int), intent(in), value :: nbr_grids
        TYPE(c_ptr), intent(out) :: grid_names(nbr_grids)
      end subroutine yac_cget_comp_grid_names_c
   end interface

   type(yac_string), allocatable :: grid_names(:)
   CHARACTER(len=*), intent(in)  :: comp_name
   integer :: nbr_grids
   INTEGER :: i
   TYPE(c_ptr), allocatable :: grid_ptr(:)

   nbr_grids = yac_cget_comp_nbr_grids_c(TRIM(comp_name) // c_null_char)
   allocate(grid_ptr(nbr_grids))
   CALL yac_cget_comp_grid_names_c(TRIM(comp_name) // c_null_char, &
                                   nbr_grids,                      &
                                   grid_ptr)
   allocate(grid_names(nbr_grids))
   DO i=1,nbr_grids
      grid_names(i)%string = yac_internal_cptr2char(grid_ptr(i))
   END DO
 end function yac_fget_comp_grid_names

 function yac_fget_comp_grid_names_instance ( yac_instance_id, comp_name) result ( grid_names )

   use mo_yac_finterface, dummy => yac_fget_comp_grid_names_instance
   use mo_yac_iso_c_helpers
   use, intrinsic :: iso_c_binding, ONLY: c_ptr, c_null_char

   implicit none

   interface
      function yac_cget_comp_nbr_grids_instance_c( yac_instance_id, &
                                                   comp_name )      &
           result( nbr_grids )                                      &
           bind(c, name='yac_cget_comp_nbr_grids_instance')
        use, intrinsic :: iso_c_binding, only: c_int, c_char
        character(kind=c_char), dimension(*), intent(in) :: comp_name
        integer(kind=c_int), value, intent(in) :: yac_instance_id
        integer(kind=c_int) :: nbr_grids
      end function yac_cget_comp_nbr_grids_instance_c

      subroutine yac_cget_comp_grid_names_instance_c( yac_instance_id, &
                                                      comp_name,       &
                                                      nbr_grids,       &
                                                      grid_names )     &
           bind(c, name='yac_cget_comp_grid_names_instance')
        use, intrinsic :: iso_c_binding, ONLY: c_int, c_ptr, c_char
        character(kind=c_char), dimension(*), intent(in) :: comp_name
        integer(kind=c_int), intent(in), value :: yac_instance_id
        integer(kind=c_int), intent(in), value :: nbr_grids
        TYPE(c_ptr), intent(out) :: grid_names(nbr_grids)
      end subroutine yac_cget_comp_grid_names_instance_c
   end interface

   integer, intent(in) :: yac_instance_id
   character(len=*), intent(in) :: comp_name
   type(yac_string), allocatable :: grid_names(:)
   integer :: nbr_grids
   INTEGER :: i
   TYPE(c_ptr), allocatable :: grid_ptr(:)

   nbr_grids =                                          &
    yac_cget_comp_nbr_grids_instance_c(yac_instance_id, &
                                       TRIM(comp_name) // c_null_char)
   allocate(grid_ptr(nbr_grids))
   CALL yac_cget_comp_grid_names_instance_c(yac_instance_id,                &
                                            TRIM(comp_name) // c_null_char, &
                                            nbr_grids,                      &
                                            grid_ptr)
   allocate(grid_names(nbr_grids))
   DO i=1,nbr_grids
      grid_names(i)%string = yac_internal_cptr2char(grid_ptr(i))
   END DO
 end function yac_fget_comp_grid_names_instance

 ! ---------------------------------------------------------------------

 function yac_fget_field_names ( comp_name, grid_name ) result( field_names )

   use mo_yac_finterface, dummy => yac_fget_field_names
   use mo_yac_iso_c_helpers
   use, intrinsic :: iso_c_binding, ONLY: c_ptr, c_null_char

   implicit none

   interface
      function yac_cget_nbr_fields_c(comp_name, grid_name) &
           result( nbr_fields )                            &
           bind(c, name='yac_cget_nbr_fields')
        use, intrinsic :: iso_c_binding, only: c_int, c_char
        character(kind=c_char), dimension(*), intent(in) :: comp_name
        character(kind=c_char), dimension(*), intent(in) :: grid_name
        integer(kind=c_int) :: nbr_fields
      end function yac_cget_nbr_fields_c

      subroutine yac_cget_field_names_c( comp_name,   &
                                         grid_name,   &
                                        nbr_fields,   &
                                        field_names ) &
           bind(c, name='yac_cget_field_names')
        use, intrinsic :: iso_c_binding, ONLY: c_int, c_ptr, c_char
        character(kind=c_char), dimension(*), intent(in) :: comp_name
        character(kind=c_char), dimension(*), intent(in) :: grid_name
        integer(kind=c_int), intent(in), value :: nbr_fields
        TYPE(c_ptr), intent(out) :: field_names(nbr_fields)
      end subroutine yac_cget_field_names_c
   end interface

   character(len=*), intent(in) :: comp_name
   character(len=*), intent(in) :: grid_name
   type(yac_string), allocatable :: field_names(:)
   integer :: nbr_fields
   INTEGER :: i
   TYPE(c_ptr), allocatable :: field_ptr(:)

   nbr_fields = yac_cget_nbr_fields_c(TRIM(comp_name)//c_null_char, &
                                      TRIM(grid_name)//c_null_char)
   allocate(field_ptr(nbr_fields))
   CALL yac_cget_field_names_c(TRIM(comp_name)//c_null_char, &
                               TRIM(grid_name)//c_null_char, &
                               nbr_fields,                   &
                               field_ptr)
   allocate(field_names(nbr_fields))
   DO i=1,nbr_fields
      field_names(i)%string = yac_internal_cptr2char(field_ptr(i))
   END DO
 end function yac_fget_field_names

 function yac_fget_field_names_instance ( yac_instance_id, &
                                          comp_name,       &
                                          grid_name )      &
      result( field_names )

   use mo_yac_finterface, dummy => yac_fget_field_names_instance
   use mo_yac_iso_c_helpers
   use, intrinsic :: iso_c_binding, ONLY: c_ptr, c_null_char

   implicit none

   interface
      function yac_cget_nbr_fields_instance_c( yac_instance_id, &
                                               comp_name,       &
                                               grid_name )      &
           result( nbr_fields ) bind(c, name='yac_cget_nbr_fields_instance')
        use, intrinsic :: iso_c_binding, only: c_int, c_char
        integer(kind=c_int), value, intent(in) :: yac_instance_id
        character(kind=c_char), dimension(*), intent(in) :: comp_name
        character(kind=c_char), dimension(*), intent(in) :: grid_name
        integer(kind=c_int) :: nbr_fields
      end function yac_cget_nbr_fields_instance_c

      subroutine yac_cget_field_names_instance_c( yac_instance_id, &
                                                  comp_name,       &
                                                  grid_name,       &
                                                  nbr_fields,      &
                                                  field_names )    &
           bind(c, name='yac_cget_field_names_instance')
        use, intrinsic :: iso_c_binding, ONLY: c_int, c_ptr, c_char
        integer(kind=c_int), intent(in), value :: yac_instance_id
        character(kind=c_char), dimension(*), intent(in) :: comp_name
        character(kind=c_char), dimension(*), intent(in) :: grid_name
        integer(kind=c_int), intent(in), value :: nbr_fields
        TYPE(c_ptr), intent(out) :: field_names(nbr_fields)
      end subroutine yac_cget_field_names_instance_c
   end interface

   integer, intent(in) :: yac_instance_id
   character(len=*), intent(in) :: comp_name
   character(len=*), intent(in) :: grid_name
   type(yac_string), allocatable :: field_names(:)
   integer :: nbr_fields
   INTEGER :: i
   TYPE(c_ptr), allocatable :: field_ptr(:)

   nbr_fields =                                                   &
     yac_cget_nbr_fields_instance_c(yac_instance_id,              &
                                    TRIM(comp_name)//c_null_char, &
                                    TRIM(grid_name)//c_null_char)
   allocate(field_ptr(nbr_fields))
   CALL yac_cget_field_names_instance_c(yac_instance_id,              &
                                        TRIM(comp_name)//c_null_char, &
                                        TRIM(grid_name)//c_null_char, &
                                        nbr_fields,                   &
                                        field_ptr)
   allocate(field_names(nbr_fields))
   DO i=1,nbr_fields
      field_names(i)%string = yac_internal_cptr2char(field_ptr(i))
   END DO
 end function yac_fget_field_names_instance

! ---------------------------------------------------------------------

function yac_fget_field_id ( comp_name, grid_name, field_name ) &
     result(field_id)

  use mo_yac_finterface, dummy => yac_fget_field_id
  use, intrinsic :: iso_c_binding, only : c_null_char

  implicit none

  interface

     function yac_cget_field_id_c ( comp_name, grid_name, field_name ) &
          result(field_id) &
          bind ( c, name='yac_cget_field_id' )

       use, intrinsic :: iso_c_binding, only : c_int, c_char

       character ( kind=c_char ), dimension(*) :: comp_name
       character ( kind=c_char ), dimension(*) :: grid_name
       character ( kind=c_char ), dimension(*) :: field_name
       integer ( kind=c_int ) :: field_id

     end function yac_cget_field_id_c

  end interface

  character(len=*), intent (in)  :: comp_name
  character(len=*), intent (in)  :: grid_name
  character(len=*), intent (in)  :: field_name
  integer :: field_id

  field_id = yac_cget_field_id_c ( TRIM(comp_name)//c_null_char, &
                                   TRIM(grid_name)//c_null_char, &
                                   TRIM(field_name)//c_null_char )

end function yac_fget_field_id

function yac_fget_field_id_instance ( yac_id,      &
                                      comp_name,   &
                                      grid_name,   &
                                      field_name ) &
     result(field_id)

  use mo_yac_finterface, dummy => yac_fget_field_id_instance
  use, intrinsic :: iso_c_binding, only : c_null_char

  implicit none

  interface

     function yac_cget_field_id_instance_c ( yac_id,      &
                                             comp_name,   &
                                             grid_name,   &
                                             field_name ) &
          result(field_id)                                &
          bind ( c, name='yac_cget_field_id_instance' )

       use, intrinsic :: iso_c_binding, only : c_int, c_char

       integer( kind=c_int ), value, intent(in) :: yac_id
       character ( kind=c_char ), dimension(*) :: comp_name
       character ( kind=c_char ), dimension(*) :: grid_name
       character ( kind=c_char ), dimension(*) :: field_name
       integer ( kind=c_int )  :: field_id

     end function yac_cget_field_id_instance_c

  end interface

  integer, intent(in) :: yac_id
  character(len=*), intent (in)  :: comp_name
  character(len=*), intent (in)  :: grid_name
  character(len=*), intent (in)  :: field_name
  integer :: field_id

  field_id =                                                     &
    yac_cget_field_id_instance_c ( yac_id,                       &
                                   TRIM(comp_name)//c_null_char, &
                                   TRIM(grid_name)//c_null_char, &
                                   TRIM(field_name)//c_null_char )

end function yac_fget_field_id_instance

! ---------------------------------------------------------------------

subroutine yac_fget_action ( field_id, action )

   use mo_yac_finterface, dummy => yac_fget_action

   implicit none

   interface

    subroutine yac_cget_action_c ( field_id, action ) &
      bind ( c, name='yac_cget_action' )

      use, intrinsic :: iso_c_binding, only : c_int

      integer ( kind=c_int ), value          :: field_id
      integer ( kind=c_int)                  :: action

    end subroutine yac_cget_action_c

   end interface

  integer, intent (in)  :: field_id !< [IN]  field identifier
  integer, intent (out) :: action   !< [OUT] action for the current timestep\n
                                    !!       (\ref YAC_ACTION_NONE,
                                    !!        \ref YAC_ACTION_COUPLING,
                                    !!        \ref YAC_ACTION_GET_FOR_RESTART,
                                    !!        \ref YAC_ACTION_PUT_FOR_RESTART, or
                                    !!        \ref YAC_ACTION_OUT_OF_BOUND)

  call yac_cget_action_c(field_id, action)

end subroutine yac_fget_action

! ---------------------------------------------------------------------

subroutine yac_fupdate ( field_id )

   use mo_yac_finterface, dummy => yac_fupdate

   implicit none

   interface

    subroutine yac_cupdate_c ( field_id ) &
      bind ( c, name='yac_cupdate' )

      use, intrinsic :: iso_c_binding, only : c_int

      integer ( kind=c_int ), value :: field_id

    end subroutine yac_cupdate_c

   end interface

  integer, intent (in)  :: field_id !< [IN]  field identifier

  call yac_cupdate_c(field_id)

end subroutine yac_fupdate

subroutine yac_fdef_couple( src_comp_name, src_grid_name, src_field_name, &
  tgt_comp_name, tgt_grid_name, tgt_field_name, &
  coupling_timestep, time_unit, time_reduction, interp_stack_config_id, &
  src_lag, tgt_lag, weight_file, mapping_side, scale_factor, scale_summand, &
  src_mask_names, tgt_mask_name)

  use, intrinsic :: iso_c_binding, only : c_null_char, c_ptr, c_null_ptr, c_loc
  use mo_yac_finterface, dummy => yac_fdef_couple
  use mo_yac_iso_c_helpers

  implicit none

  interface

     subroutine yac_cdef_couple__c( src_comp_name,          &
                                    src_grid_name,          &
                                    src_field_name,         &
                                    tgt_comp_name,          &
                                    tgt_grid_name,          &
                                    tgt_field_name,         &
                                    coupling_timestep,      &
                                    time_unit,              &
                                    time_reduction,         &
                                    interp_stack_config_id, &
                                    src_lag,                &
                                    tgt_lag,                &
                                    weight_file,            &
                                    mapping_side,           &
                                    scale_factor,           &
                                    scale_summand,          &
                                    num_src_mask_names,     &
                                    src_mask_names,         &
                                    tgt_mask_name)          &
       bind ( c, name='yac_cdef_couple_' )

       use, intrinsic :: iso_c_binding, only : c_int, c_char, c_ptr, c_double

       character ( kind=c_char ), dimension(*) :: src_comp_name
       character ( kind=c_char ), dimension(*) :: src_grid_name
       character ( kind=c_char ), dimension(*) :: src_field_name
       character ( kind=c_char ), dimension(*) :: tgt_comp_name
       character ( kind=c_char ), dimension(*) :: tgt_grid_name
       character ( kind=c_char ), dimension(*) :: tgt_field_name
       character ( kind=c_char ), dimension(*) :: coupling_timestep
       integer ( kind=c_int ), value           :: time_unit
       integer ( kind=c_int ), value           :: time_reduction
       integer ( kind=c_int ), value           :: interp_stack_config_id
       integer ( kind=c_int ), value           :: src_lag
       integer ( kind=c_int ), value           :: tgt_lag
       type ( c_ptr ), value                   :: weight_file
       integer ( kind=c_int ), value           :: mapping_side
       real ( kind=c_double ), value           :: scale_factor
       real ( kind=c_double ), value           :: scale_summand
       integer ( kind=c_int ), value           :: num_src_mask_names
       type ( c_ptr )                          :: src_mask_names(*)
       type ( c_ptr ), value                   :: tgt_mask_name
     end subroutine yac_cdef_couple__c

  end interface

  character ( len=* ), intent(in)           :: src_comp_name
  character ( len=* ), intent(in)           :: src_grid_name
  character ( len=* ), intent(in)           :: src_field_name
  character ( len=* ), intent(in)           :: tgt_comp_name
  character ( len=* ), intent(in)           :: tgt_grid_name
  character ( len=* ), intent(in)           :: tgt_field_name
  character ( len=* ), intent(in)           :: coupling_timestep
  integer, intent(in)                       :: time_unit
  integer, intent(in)                       :: time_reduction
  integer, intent(in)                       :: interp_stack_config_id
  integer, intent(in), optional             :: src_lag
  integer, intent(in), optional             :: tgt_lag
  character ( len=* ), intent(in), optional :: weight_file
  integer, intent(in), optional             :: mapping_side
  double precision, intent(in), optional    :: scale_factor
  double precision, intent(in), optional    :: scale_summand
  type(yac_string), intent(in), optional    :: src_mask_names(:)
  character ( len=* ), intent(in), optional :: tgt_mask_name

  integer :: i, j
  integer :: src_lag_cpy, tgt_lag_cpy, mapping_side_cpy
  character(kind=c_char), target :: weight_file_cpy(YAC_MAX_CHARLEN+1)
  type(c_ptr) :: weight_file_ptr
  double precision :: scale_factor_cpy, scale_summand_cpy
  integer :: num_src_mask_names
  character(kind=c_char), allocatable, target :: src_mask_names_cpy(:,:)
  type(c_ptr), allocatable :: src_mask_names_ptr(:)
  character(kind=c_char), target :: tgt_mask_name_cpy(YAC_MAX_CHARLEN+1)
  type(c_ptr) :: tgt_mask_name_ptr

  call yac_check_strlength ( src_comp_name )
  call yac_check_strlength ( src_grid_name )
  call yac_check_strlength ( src_field_name )
  call yac_check_strlength ( tgt_comp_name )
  call yac_check_strlength ( tgt_grid_name )
  call yac_check_strlength ( tgt_field_name )
  call yac_check_strlength ( coupling_timestep )
  if ( present(src_lag) ) then
    src_lag_cpy = src_lag
  else
    src_lag_cpy = 0
  end if
  if ( present(tgt_lag) ) then
    tgt_lag_cpy = tgt_lag
  else
    tgt_lag_cpy = 0
  end if
  if ( present(weight_file) ) then
     call yac_check_strlength ( weight_file )
     weight_file_cpy = c_null_char
     do i = 1, len_trim(weight_file)
       weight_file_cpy(i) = weight_file(i:i)
     end do
     weight_file_ptr = c_loc(weight_file_cpy(1))
  else
     weight_file_ptr = c_null_ptr
  end if
  if ( present(mapping_side) ) then
    mapping_side_cpy = mapping_side
  else
    mapping_side_cpy = 1
  end if
  if ( present(scale_factor) ) then
    scale_factor_cpy = scale_factor
  else
    scale_factor_cpy = 1.0
  end if
  if ( present(scale_summand) ) then
    scale_summand_cpy = scale_summand
  else
    scale_summand_cpy = 0.0
  end if
  if ( present(src_mask_names) ) then
    num_src_mask_names = size(src_mask_names)
    allocate(src_mask_names_ptr(num_src_mask_names))
    allocate(src_mask_names_cpy(YAC_MAX_CHARLEN+1,num_src_mask_names))
    src_mask_names_cpy = c_null_char
    do i = 1, num_src_mask_names
      YAC_FASSERT(allocated(src_mask_names(i)%string), "ERROR(yac_fdef_couple): source mask name not allocated")
      call yac_check_strlength ( src_mask_names(i)%string )
      do j = 1, len_trim(src_mask_names(i)%string)
        src_mask_names_cpy(j, i) = src_mask_names(i)%string(j:j)
      end do
      src_mask_names_ptr(i) = c_loc(src_mask_names_cpy(1,i))
    end do
  else
    num_src_mask_names = 0
    allocate(src_mask_names_ptr(0))
  end if
  if ( present(tgt_mask_name) ) then
     call yac_check_strlength ( tgt_mask_name )
     tgt_mask_name_cpy = c_null_char
     do i = 1, len_trim(tgt_mask_name)
       tgt_mask_name_cpy(i) = tgt_mask_name(i:i)
     end do
     tgt_mask_name_ptr = c_loc(tgt_mask_name_cpy(1))
  else
     tgt_mask_name_ptr = c_null_ptr
  end if

  call yac_cdef_couple__c( TRIM(src_comp_name) // c_null_char,     &
                           TRIM(src_grid_name) // c_null_char,     &
                           TRIM(src_field_name) // c_null_char,    &
                           TRIM(tgt_comp_name) // c_null_char,     &
                           TRIM(tgt_grid_name) // c_null_char,     &
                           TRIM(tgt_field_name) // c_null_char,    &
                           TRIM(coupling_timestep) // c_null_char, &
                           time_unit,                              &
                           time_reduction,                         &
                           interp_stack_config_id,                 &
                           src_lag_cpy,                            &
                           tgt_lag_cpy,                            &
                           weight_file_ptr,                        &
                           mapping_side_cpy,                       &
                           scale_factor_cpy,                       &
                           scale_summand_cpy,                      &
                           num_src_mask_names,                     &
                           src_mask_names_ptr,                     &
                           tgt_mask_name_ptr)

end subroutine yac_fdef_couple

subroutine yac_fdef_couple_instance( instance_id,  &
     src_comp_name, src_grid_name, src_field_name, &
     tgt_comp_name, tgt_grid_name, tgt_field_name, &
     coupling_timestep, time_unit, time_reduction, &
     interp_stack_config_id, src_lag, tgt_lag,     &
     weight_file, mapping_side, scale_factor,      &
     scale_summand, src_mask_names, tgt_mask_name )

  use, intrinsic :: iso_c_binding, only : c_null_char, c_ptr, c_null_ptr, c_loc
  use mo_yac_finterface, dummy => yac_fdef_couple_instance
  use mo_yac_iso_c_helpers

  implicit none

  interface

     subroutine yac_cdef_couple_instance__c( instance_id,            &
                                             src_comp_name,          &
                                             src_grid_name,          &
                                             src_field_name,         &
                                             tgt_comp_name,          &
                                             tgt_grid_name,          &
                                             tgt_field_name,         &
                                             coupling_timestep,      &
                                             time_unit,              &
                                             time_reduction,         &
                                             interp_stack_config_id, &
                                             src_lag,                &
                                             tgt_lag,                &
                                             weight_file,            &
                                             mapping_side,           &
                                             scale_factor,           &
                                             scale_summand,          &
                                             num_src_mask_names,     &
                                             src_mask_names,         &
                                             tgt_mask_name)          &
          bind ( c, name='yac_cdef_couple_instance_' )

       use, intrinsic :: iso_c_binding, only : c_int, c_char, c_ptr, c_double

       integer ( kind=c_int ), value           :: instance_id
       character ( kind=c_char ), dimension(*) :: src_comp_name
       character ( kind=c_char ), dimension(*) :: src_grid_name
       character ( kind=c_char ), dimension(*) :: src_field_name
       character ( kind=c_char ), dimension(*) :: tgt_comp_name
       character ( kind=c_char ), dimension(*) :: tgt_grid_name
       character ( kind=c_char ), dimension(*) :: tgt_field_name
       character ( kind=c_char ), dimension(*) :: coupling_timestep
       integer ( kind=c_int ), value           :: time_unit
       integer ( kind=c_int ), value           :: time_reduction
       integer ( kind=c_int ), value           :: interp_stack_config_id
       integer ( kind=c_int ), value           :: src_lag
       integer ( kind=c_int ), value           :: tgt_lag
       type ( c_ptr ), value                   :: weight_file
       integer ( kind=c_int ), value           :: mapping_side
       real ( kind=c_double ), value           :: scale_factor
       real ( kind=c_double ), value           :: scale_summand
       integer ( kind=c_int ), value           :: num_src_mask_names
       type ( c_ptr )                          :: src_mask_names(*)
       type ( c_ptr ), value                   :: tgt_mask_name
     end subroutine yac_cdef_couple_instance__c

  end interface

  integer, intent(in)                       :: instance_id
  character ( len=* ), intent(in)           :: src_comp_name
  character ( len=* ), intent(in)           :: src_grid_name
  character ( len=* ), intent(in)           :: src_field_name
  character ( len=* ), intent(in)           :: tgt_comp_name
  character ( len=* ), intent(in)           :: tgt_grid_name
  character ( len=* ), intent(in)           :: tgt_field_name
  character ( len=* ), intent(in)           :: coupling_timestep
  integer, intent(in)                       :: time_unit
  integer, intent(in)                       :: time_reduction
  integer, intent(in)                       :: interp_stack_config_id
  integer, intent(in), optional             :: src_lag
  integer, intent(in), optional             :: tgt_lag
  character ( len=* ), intent(in), optional :: weight_file
  integer, intent(in), optional             :: mapping_side
  double precision, intent(in), optional    :: scale_factor
  double precision, intent(in), optional    :: scale_summand
  type(yac_string), intent(in), optional    :: src_mask_names(:)
  character ( len=* ), intent(in), optional :: tgt_mask_name

  integer :: i, j
  integer :: src_lag_cpy, tgt_lag_cpy, mapping_side_cpy
  character(kind=c_char), target :: weight_file_cpy(YAC_MAX_CHARLEN+1)
  type(c_ptr) :: weight_file_ptr
  double precision :: scale_factor_cpy, scale_summand_cpy
  integer :: num_src_mask_names
  character(kind=c_char), allocatable, target :: src_mask_names_cpy(:,:)
  type(c_ptr), allocatable :: src_mask_names_ptr(:)
  character(kind=c_char), target :: tgt_mask_name_cpy(YAC_MAX_CHARLEN+1)
  type(c_ptr) :: tgt_mask_name_ptr

  call yac_check_strlength ( src_comp_name )
  call yac_check_strlength ( src_grid_name )
  call yac_check_strlength ( src_field_name )
  call yac_check_strlength ( tgt_comp_name )
  call yac_check_strlength ( tgt_grid_name )
  call yac_check_strlength ( tgt_field_name )
  call yac_check_strlength ( coupling_timestep )
  if ( present(src_lag) ) then
    src_lag_cpy = src_lag
  else
    src_lag_cpy = 0
  end if
  if ( present(tgt_lag) ) then
    tgt_lag_cpy = tgt_lag
  else
    tgt_lag_cpy = 0
  end if
  if ( present(weight_file) ) then
     call yac_check_strlength ( weight_file )
     weight_file_cpy = c_null_char
     do i = 1, len_trim(weight_file)
       weight_file_cpy(i) = weight_file(i:i)
     end do
     weight_file_ptr = c_loc(weight_file_cpy(1))
  else
     weight_file_ptr = c_null_ptr
  end if
  if ( present(mapping_side) ) then
    mapping_side_cpy = mapping_side
  else
    mapping_side_cpy = 1
  end if
  if ( present(scale_factor) ) then
    scale_factor_cpy = scale_factor
  else
    scale_factor_cpy = 1.0
  end if
  if ( present(scale_summand) ) then
    scale_summand_cpy = scale_summand
  else
    scale_summand_cpy = 0.0
  end if
  if ( present(src_mask_names) ) then
    num_src_mask_names = size(src_mask_names)
    allocate(src_mask_names_ptr(num_src_mask_names))
    allocate(src_mask_names_cpy(YAC_MAX_CHARLEN+1,num_src_mask_names))
    src_mask_names_cpy = c_null_char
    do i = 1, num_src_mask_names
      YAC_FASSERT(allocated(src_mask_names(i)%string), "ERROR(yac_fdef_couple): source mask name not allocated")
      call yac_check_strlength ( src_mask_names(i)%string )
      do j = 1, len_trim(src_mask_names(i)%string)
        src_mask_names_cpy(j, i) = src_mask_names(i)%string(j:j)
      end do
      src_mask_names_ptr(i) = c_loc(src_mask_names_cpy(1,i))
    end do
  else
    num_src_mask_names = 0
    allocate(src_mask_names_ptr(0))
  end if
  if ( present(tgt_mask_name) ) then
     call yac_check_strlength ( tgt_mask_name )
     tgt_mask_name_cpy = c_null_char
     do i = 1, len_trim(tgt_mask_name)
       tgt_mask_name_cpy(i) = tgt_mask_name(i:i)
     end do
     tgt_mask_name_ptr = c_loc(tgt_mask_name_cpy(1))
  else
     tgt_mask_name_ptr = c_null_ptr
  end if

  call yac_cdef_couple_instance__c( instance_id, &
                                    TRIM(src_comp_name) // c_null_char,     &
                                    TRIM(src_grid_name) // c_null_char,     &
                                    TRIM(src_field_name) // c_null_char,    &
                                    TRIM(tgt_comp_name) // c_null_char,     &
                                    TRIM(tgt_grid_name) // c_null_char,     &
                                    TRIM(tgt_field_name) // c_null_char,    &
                                    TRIM(coupling_timestep) // c_null_char, &
                                    time_unit,                              &
                                    time_reduction,                         &
                                    interp_stack_config_id,                 &
                                    src_lag_cpy,                            &
                                    tgt_lag_cpy,                            &
                                    weight_file_ptr,                        &
                                    mapping_side_cpy,                       &
                                    scale_factor_cpy,                       &
                                    scale_summand_cpy,                      &
                                    num_src_mask_names,                     &
                                    src_mask_names_ptr,                     &
                                    tgt_mask_name_ptr )

end subroutine yac_fdef_couple_instance

! ---------------------------------------------------------------------

function yac_fget_component_name_from_field_id ( field_id ) &
     result(comp_name)

   use mo_yac_finterface, dummy => yac_fget_component_name_from_field_id
   use mo_yac_iso_c_helpers
   use, intrinsic :: iso_c_binding, only: c_null_char

   implicit none

   interface

      function yac_cget_component_name_from_field_id_c ( field_id ) &
           result(comp_name)                                        &
       bind ( c, name='yac_cget_component_name_from_field_id' )

       use, intrinsic :: iso_c_binding, only : c_int, c_ptr

       integer ( kind=c_int ), value :: field_id !< [IN]  field ID
       TYPE(c_ptr) :: comp_name

     end function yac_cget_component_name_from_field_id_c

   end interface

   integer, intent (in)           :: field_id  !< [IN]  field identifier
   character (len=:), allocatable :: comp_name !< [OUT] field name

   comp_name = yac_internal_cptr2char( &
        yac_cget_component_name_from_field_id_c ( field_id ))

 end function yac_fget_component_name_from_field_id

 ! ---------------------------------------------------------------------

function yac_fget_grid_name_from_field_id ( field_id ) &
     result(grid_name)

   use mo_yac_finterface, dummy => yac_fget_grid_name_from_field_id
   use mo_yac_iso_c_helpers
   use, intrinsic :: iso_c_binding, only: c_null_char

   implicit none

   interface

      function yac_cget_grid_name_from_field_id_c ( field_id ) &
           result(grid_name)                                   &
       bind ( c, name='yac_cget_grid_name_from_field_id' )

       use, intrinsic :: iso_c_binding, only : c_int, c_ptr

       integer ( kind=c_int ), value  :: field_id !< [IN]  field ID
       TYPE(c_ptr) :: grid_name

     end function yac_cget_grid_name_from_field_id_c

   end interface

   integer, intent (in)           :: field_id  !< [IN]  field identifier
   character (len=:), ALLOCATABLE :: grid_name !< [OUT] field name

   grid_name = yac_internal_cptr2char( &
        yac_cget_grid_name_from_field_id_c ( field_id ))

 end function yac_fget_grid_name_from_field_id

! ---------------------------------------------------------------------

function yac_fget_field_name_from_field_id ( field_id ) &
     result(field_name)

   use mo_yac_finterface, dummy => yac_fget_field_name_from_field_id
   use mo_yac_iso_c_helpers
   use, intrinsic :: iso_c_binding, only: c_null_char

   implicit none

   interface

      function yac_cget_field_name_from_field_id_c ( field_id ) &
           result(field_name)                                   &
       bind ( c, name='yac_cget_field_name_from_field_id' )

       use, intrinsic :: iso_c_binding, only : c_int, c_ptr

       integer ( kind=c_int ), value :: field_id !< [IN]  field ID
       TYPE(c_ptr)                   :: field_name

     end function yac_cget_field_name_from_field_id_c

   end interface

   integer, intent (in)           :: field_id   !< [IN]  field identifier
   character (len=:), ALLOCATABLE :: field_name !< [OUT] field name

   field_name = yac_internal_cptr2char( &
        yac_cget_field_name_from_field_id_c ( field_id ))

 end function yac_fget_field_name_from_field_id

! ---------------------------------------------------------------------

 function yac_fget_role_from_field_id ( field_id )

   use mo_yac_finterface, dummy => yac_fget_role_from_field_id
   use, intrinsic :: iso_c_binding, only: c_null_char

   implicit none

   interface

     function yac_cget_role_from_field_id_c ( field_id ) &
       bind ( c, name='yac_cget_role_from_field_id' )

       use, intrinsic :: iso_c_binding, only : c_int, c_char

       integer ( kind=c_int ), value :: field_id   !< [IN]  field ID
       integer ( kind=c_int )        :: yac_cget_role_from_field_id_c

     end function  yac_cget_role_from_field_id_c

   end interface

   integer, intent (in) :: field_id   !< [IN]  field identifier
   integer              :: yac_fget_role_from_field_id

   yac_fget_role_from_field_id = &
     yac_cget_role_from_field_id_c ( field_id )

 end function yac_fget_role_from_field_id

 ! Note that in contrast to most of the other functions in this file, we have to
 ! introduce a separate result variable below. Otherwise, NVHPC fails to compile
 ! the file due to the fact that the name of the function is also a name of the
 ! interface in mo_yac_finterface.
 function yac_fget_field_role ( comp_name, grid_name, field_name ) result( res )

   use mo_yac_finterface, dummy => yac_fget_field_role
   use, intrinsic :: iso_c_binding, only: c_null_char

   implicit none

   interface

     function yac_cget_field_role_c ( comp_name, grid_name, field_name ) &
       bind ( c, name='yac_cget_field_role' )

       use, intrinsic :: iso_c_binding, only : c_int, c_char

       character( kind=c_char), dimension(*), intent(in)  :: comp_name
       character( kind=c_char), dimension(*), intent(in)  :: grid_name
       character( kind=c_char), dimension(*), intent(in)  :: field_name
       integer ( kind=c_int ) :: yac_cget_field_role_c

     end function  yac_cget_field_role_c

   end interface

   character(len=*), intent(in) :: comp_name
   character(len=*), intent(in) :: grid_name
   character(len=*), intent(in) :: field_name
   integer                      :: res

   res = yac_cget_field_role_c (     &
     TRIM(comp_name) // c_null_char, &
     TRIM(grid_name) // c_null_char, &
     TRIM(field_name) // c_null_char )

 end function yac_fget_field_role

 function yac_fget_field_role_instance ( yac_instance_id, comp_name, grid_name, field_name )

   use mo_yac_finterface, dummy => yac_fget_field_role_instance
   use, intrinsic :: iso_c_binding, only: c_null_char

   implicit none

   interface

     function yac_cget_field_role_instance_c ( yac_instance_id, &
                                               comp_name,       &
                                               grid_name,       &
                                               field_name )     &
       bind ( c, name='yac_cget_field_role_instance' )

       use, intrinsic :: iso_c_binding, only : c_int, c_char

       integer( kind=c_int ), intent(in), value :: yac_instance_id
       character( kind=c_char), dimension(*), intent(in)  :: comp_name
       character( kind=c_char), dimension(*), intent(in)  :: grid_name
       character( kind=c_char), dimension(*), intent(in)  :: field_name
       integer ( kind=c_int ) :: yac_cget_field_role_instance_c

     end function  yac_cget_field_role_instance_c

   end interface

   integer, intent(in) :: yac_instance_id
   character(len=*), intent(in) :: comp_name
   character(len=*), intent(in) :: grid_name
   character(len=*), intent(in) :: field_name
   integer :: yac_fget_field_role_instance

   yac_fget_field_role_instance =       &
     yac_cget_field_role_instance_c (   &
        yac_instance_id,                &
        TRIM(comp_name) // c_null_char, &
        TRIM(grid_name) // c_null_char, &
        TRIM(field_name) // c_null_char )

 end function yac_fget_field_role_instance

! ---------------------------------------------------------------------

function yac_fget_timestep_from_field_id ( field_id ) result(string)

   use mo_yac_finterface, dummy => yac_fget_timestep_from_field_id
   use mo_yac_iso_c_helpers

   implicit none

   interface

      function yac_cget_timestep_from_field_id_c ( field_id ) &
           result(string)                                     &
           bind ( c, name='yac_cget_timestep_from_field_id' )

       use, intrinsic :: iso_c_binding, only : c_int, c_ptr

       integer ( kind=c_int ), value :: field_id
       type(c_ptr)                   :: string

     end function  yac_cget_timestep_from_field_id_c

    end interface

    integer, intent (in)           :: field_id !< [IN]  field identifier
    character (len=:), ALLOCATABLE :: string   !< [OUT] timestep in iso format

    string =                  &
      yac_internal_cptr2char( &
        yac_cget_timestep_from_field_id_c(field_id))

  end function yac_fget_timestep_from_field_id

  function yac_fget_field_timestep ( comp_name, grid_name, field_name ) &
       result( timestep )

    use mo_yac_finterface, dummy => yac_fget_field_timestep
    use mo_yac_iso_c_helpers
    use, intrinsic :: iso_c_binding, only: c_ptr, c_null_char

   implicit none

   interface

      function yac_cget_field_timestep_c ( comp_name, grid_name, field_name) &
           result( timestep ) &
           bind ( c, name='yac_cget_field_timestep' )

        use, intrinsic :: iso_c_binding, only : c_int, c_char, c_ptr

       character( kind=c_char), dimension(*), intent(in)  :: comp_name
       character( kind=c_char), dimension(*), intent(in)  :: grid_name
       character( kind=c_char), dimension(*), intent(in)  :: field_name
       type(c_ptr)  :: timestep
     end function  yac_cget_field_timestep_c

   end interface

   character(len=*), intent(in)  :: comp_name
   character(len=*), intent(in)  :: grid_name
   character(len=*), intent(in)  :: field_name
   character(len=:), ALLOCATABLE :: timestep
   TYPE(c_ptr)                   :: c_char_ptr

   c_char_ptr =                         &
     yac_cget_field_timestep_c (        &
        TRIM(comp_name) // c_null_char, &
        TRIM(grid_name) // c_null_char, &
        TRIM(field_name) // c_null_char)

   timestep = yac_internal_cptr2char(c_char_ptr)
 end function  yac_fget_field_timestep

 function yac_fget_field_timestep_instance ( yac_instance_id, comp_name, grid_name, field_name ) &
      result( timestep )

   use mo_yac_finterface, dummy => yac_fget_field_timestep_instance
   use mo_yac_iso_c_helpers
   use, intrinsic :: iso_c_binding, only: c_ptr, c_null_char

   implicit none

   interface

      function yac_cget_field_timestep_instance_c ( yac_instance_id, &
                                                    comp_name,       &
                                                    grid_name,       &
                                                    field_name )     &
           result( timestep )                                        &
           bind ( c, name='yac_cget_field_timestep_instance' )

        use, intrinsic :: iso_c_binding, only : c_int, c_char, c_ptr

        integer( kind=c_int ), intent(in), value :: yac_instance_id
        character( kind=c_char), dimension(*), intent(in) :: comp_name
        character( kind=c_char), dimension(*), intent(in) :: grid_name
        character( kind=c_char), dimension(*), intent(in) :: field_name
        type(c_ptr)  :: timestep

     end function  yac_cget_field_timestep_instance_c

   end interface

   integer, intent(in) :: yac_instance_id
   character(len=*), intent(in) :: comp_name
   character(len=*), intent(in) :: grid_name
   character(len=*), intent(in) :: field_name
   character(len=:), ALLOCATABLE :: timestep
   type(c_ptr)                  :: c_char_ptr

   c_char_ptr =                           &
     yac_cget_field_timestep_instance_c ( &
        yac_instance_id,                  &
        TRIM(comp_name) // c_null_char,   &
        TRIM(grid_name) // c_null_char,   &
        TRIM(field_name) // c_null_char)

   timestep = yac_internal_cptr2char(c_char_ptr)

 end function   yac_fget_field_timestep_instance

! ---------------------------------------------------------------------

function yac_fget_field_datetime(field_id) &
     result(datetime)

  use mo_yac_finterface, dummy => yac_fget_field_datetime
  use mo_yac_iso_c_helpers
  use, intrinsic :: iso_c_binding, only : c_char, c_ptr

  implicit none

  interface
     function yac_cget_field_datetime_c( field_id ) &
          result(datetime) &
          bind ( c, name="yac_cget_field_datetime" )
       use, intrinsic :: iso_c_binding, only : c_int, c_char, c_ptr
       integer (kind=c_int ), intent(in), value :: field_id
       type(c_ptr) :: datetime
     end function yac_cget_field_datetime_c
  end interface

  integer, intent(in) :: field_id
  character(len=:), allocatable :: datetime
  datetime = yac_internal_cptr2char(yac_cget_field_datetime_c(field_id))
end function yac_fget_field_datetime

! ---------------------------------------------------------------------

 subroutine yac_fenable_field_frac_mask( &
   comp_name, grid_name, field_name, frac_mask_fallback_value)

   use, intrinsic :: iso_c_binding, only: c_null_char
   use mo_yac_finterface, dummy => yac_fenable_field_frac_mask

   implicit none

   interface
      subroutine yac_cenable_field_frac_mask_c (comp_name,                &
                                                grid_name,                &
                                                field_name,               &
                                                frac_mask_fallback_value) &
           bind( c, name="yac_cenable_field_frac_mask" )
        use, intrinsic :: iso_c_binding, only: c_char, c_double
        character(kind=c_char), dimension(*), intent(in) :: comp_name
        character(kind=c_char), dimension(*), intent(in) :: grid_name
        character(kind=c_char), dimension(*), intent(in) :: field_name
        real(kind=c_double), value :: frac_mask_fallback_value
      end subroutine yac_cenable_field_frac_mask_c
   end interface

   character(len=*), intent(in) :: comp_name
   character(len=*), intent(in) :: grid_name
   character(len=*), intent(in) :: field_name
   double precision, intent(in) :: frac_mask_fallback_value

   CALL yac_cenable_field_frac_mask_c (  &
        TRIM(comp_name) // c_null_char,  &
        TRIM(grid_name) // c_null_char,  &
        TRIM(field_name) // c_null_char, &
        frac_mask_fallback_value)
 end subroutine yac_fenable_field_frac_mask

 subroutine yac_fenable_field_frac_mask_instance(yac_instance_id, comp_name, &
      grid_name, field_name, frac_mask_fallback_value)

   use, intrinsic :: iso_c_binding, only: c_null_char
   use mo_yac_finterface, dummy => yac_fenable_field_frac_mask_instance

   implicit none

   interface
      subroutine yac_cenable_field_frac_mask_instance_c ( &
        yac_instance_id, comp_name, grid_name, field_name, &
        frac_mask_fallback_value) &
          bind( c, name="yac_cenable_field_frac_mask_instance" )
        use, intrinsic :: iso_c_binding, only: c_char, c_double, c_int
        integer(kind=c_int), intent(in), value           :: yac_instance_id
        character(kind=c_char), dimension(*), intent(in) :: comp_name
        character(kind=c_char), dimension(*), intent(in) :: grid_name
        character(kind=c_char), dimension(*), intent(in) :: field_name
        real(kind=c_double), value :: frac_mask_fallback_value
      end subroutine yac_cenable_field_frac_mask_instance_c
   end interface

   integer, intent(in)          :: yac_instance_id
   character(len=*), intent(in) :: comp_name
   character(len=*), intent(in) :: grid_name
   character(len=*), intent(in) :: field_name
   double precision, intent(in) :: frac_mask_fallback_value

   CALL yac_cenable_field_frac_mask_instance_c (  &
        yac_instance_id,                          &
        TRIM(comp_name) // c_null_char,           &
        TRIM(grid_name) // c_null_char,           &
        TRIM(field_name) // c_null_char,          &
        frac_mask_fallback_value)
 end subroutine yac_fenable_field_frac_mask_instance

! ---------------------------------------------------------------------

 ! Note that in contrast to most of the other functions in this file, we have to
 ! introduce a separate result variable below. Otherwise, NVHPC fails to compile
 ! the file due to the fact that the name of the function is also a name of the
 ! interface in mo_yac_finterface.
 function yac_fget_field_frac_mask_fallback_value ( &
   comp_name, grid_name, field_name ) result( res )

   use mo_yac_finterface, dummy => yac_fget_field_frac_mask_fallback_value
   use, intrinsic :: iso_c_binding, only: c_null_char

   implicit none

   interface

     function yac_cget_field_frac_mask_fallback_value_c ( &
       comp_name, grid_name, field_name ) &
       bind ( c, name='yac_cget_field_frac_mask_fallback_value' )

       use, intrinsic :: iso_c_binding, only : c_double, c_char

       character( kind=c_char), dimension(*), intent(in)  :: comp_name
       character( kind=c_char), dimension(*), intent(in)  :: grid_name
       character( kind=c_char), dimension(*), intent(in)  :: field_name
       real ( kind=c_double ) :: yac_cget_field_frac_mask_fallback_value_c

     end function  yac_cget_field_frac_mask_fallback_value_c

   end interface

   character(len=*), intent(in) :: comp_name
   character(len=*), intent(in) :: grid_name
   character(len=*), intent(in) :: field_name
   double precision             :: res

   res = yac_cget_field_frac_mask_fallback_value_c ( &
     TRIM(comp_name) // c_null_char,                 &
     TRIM(grid_name) // c_null_char,                 &
     TRIM(field_name) // c_null_char )

 end function yac_fget_field_frac_mask_fallback_value

 function yac_fget_field_frac_mask_fallback_value_instance ( &
   yac_instance_id, comp_name, grid_name, field_name )

   use mo_yac_finterface, dummy => yac_fget_field_frac_mask_fallback_value_instance
   use, intrinsic :: iso_c_binding, only: c_null_char

   implicit none

   interface

     function yac_cget_field_frac_mask_fallback_value_instance_c ( &
       yac_instance_id, comp_name, grid_name, field_name ) &
       bind ( c, name='yac_cget_field_frac_mask_fallback_value_instance' )

       use, intrinsic :: iso_c_binding, only : c_int, c_char, c_double

       integer( kind=c_int ), intent(in), value :: yac_instance_id
       character( kind=c_char), dimension(*), intent(in)  :: comp_name
       character( kind=c_char), dimension(*), intent(in)  :: grid_name
       character( kind=c_char), dimension(*), intent(in)  :: field_name
       real ( kind=c_double ) :: yac_cget_field_frac_mask_fallback_value_instance_c

     end function  yac_cget_field_frac_mask_fallback_value_instance_c

   end interface

   integer, intent(in) :: yac_instance_id
   character(len=*), intent(in) :: comp_name
   character(len=*), intent(in) :: grid_name
   character(len=*), intent(in) :: field_name
   double precision :: yac_fget_field_frac_mask_fallback_value_instance

   yac_fget_field_frac_mask_fallback_value_instance =     &
     yac_cget_field_frac_mask_fallback_value_instance_c ( &
        yac_instance_id,                                  &
        TRIM(comp_name) // c_null_char,                   &
        TRIM(grid_name) // c_null_char,                   &
        TRIM(field_name) // c_null_char )

 end function yac_fget_field_frac_mask_fallback_value_instance

 ! ---------------------------------------------------------------------

function yac_fget_collection_size_from_field_id ( field_id ) &
     RESULT(collection_size)

   use mo_yac_finterface, dummy => yac_fget_collection_size_from_field_id
   use mo_yac_iso_c_helpers

   implicit none

   interface

      function yac_cget_collection_size_from_field_id_c ( field_id ) &
           RESULT(collection_size) &
           bind ( c, name='yac_cget_collection_size_from_field_id' )

       use, intrinsic :: iso_c_binding, only : c_int

       integer ( kind=c_int ), value :: field_id
       integer ( kind=c_int )        :: collection_size

     end function  yac_cget_collection_size_from_field_id_c

   end interface

   integer, intent (in)  :: field_id        !< [IN]  field identifier
   integer               :: collection_size !< [OUT] collection size

   collection_size = &
     yac_cget_collection_size_from_field_id_c ( field_id )

 end function  yac_fget_collection_size_from_field_id

 ! Note that in contrast to most of the other functions in this file, we have to
 ! introduce a separate result variable below. Otherwise, NVHPC fails to compile
 ! the file due to the fact that the name of the function is also a name of the
 ! interface in mo_yac_finterface.
 function yac_fget_field_collection_size ( comp_name, grid_name, field_name ) result( res )

   use mo_yac_finterface, dummy => yac_fget_field_collection_size
   use, intrinsic :: iso_c_binding, only: c_null_char

   implicit none

   interface

     function yac_cget_field_collection_size_c ( comp_name,  &
                                                grid_name,   &
                                                field_name ) &
       bind ( c, name='yac_cget_field_collection_size' )

       use, intrinsic :: iso_c_binding, only : c_int, c_char

       character( kind=c_char), dimension(*), intent(in)  :: comp_name
       character( kind=c_char), dimension(*), intent(in)  :: grid_name
       character( kind=c_char), dimension(*), intent(in)  :: field_name
       integer ( kind=c_int ) :: yac_cget_field_collection_size_c

     end function  yac_cget_field_collection_size_c

   end interface

   character(len=*), intent(in) :: comp_name
   character(len=*), intent(in) :: grid_name
   character(len=*), intent(in) :: field_name
   integer                      :: res

   res = yac_cget_field_collection_size_c ( &
     TRIM(comp_name) // c_null_char,        &
     TRIM(grid_name) // c_null_char,        &
     TRIM(field_name) // c_null_char )

 end function yac_fget_field_collection_size

 function yac_fget_field_collection_size_instance ( yac_instance_id, comp_name, grid_name, field_name )

   use mo_yac_finterface, dummy => yac_fget_field_collection_size_instance
   use, intrinsic :: iso_c_binding, only: c_null_char

   implicit none

   interface

     function yac_cget_field_collection_size_instance_c ( yac_instance_id, &
                                                          comp_name,       &
                                                          grid_name,       &
                                                          field_name )     &
       bind ( c, name='yac_cget_field_collection_size_instance' )

       use, intrinsic :: iso_c_binding, only : c_int, c_char

       integer( kind=c_int ), intent(in), value :: yac_instance_id
       character( kind=c_char), dimension(*), intent(in)  :: comp_name
       character( kind=c_char), dimension(*), intent(in)  :: grid_name
       character( kind=c_char), dimension(*), intent(in)  :: field_name
       integer ( kind=c_int ) :: yac_cget_field_collection_size_instance_c

     end function  yac_cget_field_collection_size_instance_c

   end interface

   integer, intent(in) :: yac_instance_id
   character(len=*), intent(in) :: comp_name
   character(len=*), intent(in) :: grid_name
   character(len=*), intent(in) :: field_name
   integer :: yac_fget_field_collection_size_instance

   yac_fget_field_collection_size_instance =     &
     yac_cget_field_collection_size_instance_c ( &
        yac_instance_id,                         &
        TRIM(comp_name) // c_null_char,          &
        TRIM(grid_name) // c_null_char,          &
        TRIM(field_name) // c_null_char )

 end function yac_fget_field_collection_size_instance

 ! ---------------------------------------------------------------------

 subroutine yac_fdef_component_metadata(comp_name, metadata)

   use, intrinsic :: iso_c_binding, only: c_null_char
   use mo_yac_finterface, dummy => yac_fdef_component_metadata

   implicit none

   interface
      subroutine yac_cdef_component_metadata_c (comp_name, metadata) &
           bind( c, name="yac_cdef_component_metadata" )
        use, intrinsic :: iso_c_binding, only: c_char
        character(kind=c_char), dimension(*), intent(in) :: comp_name
        character(kind=c_char), dimension(*), intent(in) :: metadata
      end subroutine yac_cdef_component_metadata_c
   end interface

   character(len=*), intent(in) :: comp_name
   character(len=*), intent(in) :: metadata
   CALL yac_cdef_component_metadata_c (TRIM(comp_name) // c_null_char, &
                                       TRIM(metadata) // c_null_char)
 end subroutine yac_fdef_component_metadata

 subroutine yac_fdef_component_metadata_instance(yac_instance_id, comp_name, &
      metadata)

   use, intrinsic :: iso_c_binding, only: c_null_char
   use mo_yac_finterface, dummy => yac_fdef_component_metadata_instance

   implicit none

   interface
      subroutine yac_cdef_component_metadata_instance_c (yac_instance_id, &
                                                         comp_name,       &
                                                         metadata)        &
           bind( c, name="yac_cdef_component_metadata_instance" )
        use, intrinsic :: iso_c_binding, only: c_char, c_int
        integer(kind=c_int), intent(in), value           :: yac_instance_id
        character(kind=c_char), dimension(*), intent(in) :: comp_name
        character(kind=c_char), dimension(*), intent(in) :: metadata
      end subroutine yac_cdef_component_metadata_instance_c
   end interface

   integer, intent(in)          :: yac_instance_id
   character(len=*), intent(in) :: comp_name
   character(len=*), intent(in) :: metadata
   CALL yac_cdef_component_metadata_instance_c ( &
        yac_instance_id,                         &
        TRIM(comp_name) // c_null_char,          &
        TRIM(metadata) // c_null_char)
 end subroutine yac_fdef_component_metadata_instance

 subroutine yac_fdef_grid_metadata(grid_name, metadata)

   use, intrinsic :: iso_c_binding, only: c_null_char
   use mo_yac_finterface, dummy => yac_fdef_grid_metadata

   implicit none

   interface
      subroutine yac_cdef_grid_metadata_c (grid_name, metadata) &
           bind( c, name="yac_cdef_grid_metadata" )
        use, intrinsic :: iso_c_binding, only: c_char
        character(kind=c_char), dimension(*), intent(in) :: grid_name
        character(kind=c_char), dimension(*), intent(in) :: metadata
      end subroutine yac_cdef_grid_metadata_c
   end interface

   character(len=*), intent(in) :: grid_name
   character(len=*), intent(in) :: metadata
   CALL yac_cdef_grid_metadata_c ( TRIM(grid_name) // c_null_char, &
                                   TRIM(metadata) // c_null_char)
 end subroutine yac_fdef_grid_metadata

 subroutine yac_fdef_grid_metadata_instance(yac_instance_id, grid_name, metadata)

   use, intrinsic :: iso_c_binding, only: c_null_char
   use mo_yac_finterface, dummy => yac_fdef_grid_metadata_instance

   implicit none

   interface
      subroutine yac_cdef_grid_metadata_instance_c (yac_instance_id, &
                                                    grid_name,       &
                                                    metadata)        &
           bind( c, name="yac_cdef_grid_metadata_instance" )
        use, intrinsic :: iso_c_binding, only: c_char, c_int
        integer(kind=c_int), intent(in), value           :: yac_instance_id
        character(kind=c_char), dimension(*), intent(in) :: grid_name
        character(kind=c_char), dimension(*), intent(in) :: metadata
      end subroutine yac_cdef_grid_metadata_instance_c
   end interface

   integer, intent(in)          :: yac_instance_id
   character(len=*), intent(in) :: grid_name
   character(len=*), intent(in) :: metadata
   CALL yac_cdef_grid_metadata_instance_c ( &
        yac_instance_id,                    &
        TRIM(grid_name) // c_null_char,     &
        TRIM(metadata) // c_null_char)
 end subroutine yac_fdef_grid_metadata_instance

 subroutine yac_fdef_field_metadata(comp_name, grid_name, field_name, metadata)

   use, intrinsic :: iso_c_binding, only: c_null_char
   use mo_yac_finterface, dummy => yac_fdef_field_metadata

   implicit none

   interface
      subroutine yac_cdef_field_metadata_c (comp_name,  &
                                            grid_name,  &
                                            field_name, &
                                            metadata)   &
           bind( c, name="yac_cdef_field_metadata" )
        use, intrinsic :: iso_c_binding, only: c_char
        character(kind=c_char), dimension(*), intent(in) :: comp_name
        character(kind=c_char), dimension(*), intent(in) :: grid_name
        character(kind=c_char), dimension(*), intent(in) :: field_name
        character(kind=c_char), dimension(*), intent(in) :: metadata
      end subroutine yac_cdef_field_metadata_c
   end interface

   character(len=*), intent(in) :: comp_name
   character(len=*), intent(in) :: grid_name
   character(len=*), intent(in) :: field_name
   character(len=*), intent(in) :: metadata
   CALL yac_cdef_field_metadata_c (      &
        TRIM(comp_name) // c_null_char,  &
        TRIM(grid_name) // c_null_char,  &
        TRIM(field_name) // c_null_char, &
        TRIM(metadata) // c_null_char)
 end subroutine yac_fdef_field_metadata

 subroutine yac_fdef_field_metadata_instance(yac_instance_id, comp_name, &
      grid_name, field_name, metadata)

   use, intrinsic :: iso_c_binding, only: c_null_char
   use mo_yac_finterface, dummy => yac_fdef_field_metadata_instance

   implicit none

   interface
      subroutine yac_cdef_field_metadata_instance_c (yac_instance_id, &
                                                     comp_name,       &
                                                     grid_name,       &
                                                     field_name,      &
                                                     metadata)        &
           bind( c, name="yac_cdef_field_metadata_instance" )
        use, intrinsic :: iso_c_binding, only: c_char, c_int
        integer(kind=c_int), intent(in), value           :: yac_instance_id
        character(kind=c_char), dimension(*), intent(in) :: comp_name
        character(kind=c_char), dimension(*), intent(in) :: grid_name
        character(kind=c_char), dimension(*), intent(in) :: field_name
        character(kind=c_char), dimension(*), intent(in) :: metadata
      end subroutine yac_cdef_field_metadata_instance_c
   end interface

   integer, intent(in)          :: yac_instance_id
   character(len=*), intent(in) :: comp_name
   character(len=*), intent(in) :: grid_name
   character(len=*), intent(in) :: field_name
   character(len=*), intent(in) :: metadata
   CALL yac_cdef_field_metadata_instance_c ( &
        yac_instance_id,                     &
        TRIM(comp_name) // c_null_char,      &
        TRIM(grid_name) // c_null_char,      &
        TRIM(field_name) // c_null_char,     &
        TRIM(metadata) // c_null_char)
 end subroutine yac_fdef_field_metadata_instance

 function yac_fcomponent_has_metadata(comp_name) result( has_metadata )
   use, intrinsic :: iso_c_binding, only: c_null_char, c_associated
   use mo_yac_finterface, dummy => yac_fcomponent_has_metadata
   use mo_yac_iso_c_helpers
   implicit none

   interface
      function yac_cget_component_metadata_c (comp_name) &
           result(metadata)                              &
           bind( c, name="yac_cget_component_metadata")
        use, intrinsic :: iso_c_binding, only: c_char, c_ptr
        character(kind=c_char), dimension(*), intent(in) :: comp_name
        type(c_ptr) :: metadata
      end function yac_cget_component_metadata_c
   end interface
   character(len=*), intent(in) :: comp_name
   logical :: has_metadata
   has_metadata = &
    c_associated( &
      yac_cget_component_metadata_c(TRIM(comp_name) // c_null_char))
 end function yac_fcomponent_has_metadata

 function yac_fcomponent_has_metadata_instance(yac_instance_id, comp_name) &
      result( has_metadata )
   use, intrinsic :: iso_c_binding, only: c_null_char, c_associated
   use mo_yac_finterface, dummy => yac_fcomponent_has_metadata_instance
   use mo_yac_iso_c_helpers
   implicit none

   interface
      function yac_cget_component_metadata_instance_c(yac_instance_id, &
                                                      comp_name)       &
           result(metadata)                                            &
           bind( c, name="yac_cget_component_metadata_instance")
        use, intrinsic :: iso_c_binding, only: c_char, c_ptr, c_int
        integer(kind=c_int), intent(in), value :: yac_instance_id
        character(kind=c_char), dimension(*), intent(in) :: comp_name
        type(c_ptr) :: metadata
      end function yac_cget_component_metadata_instance_c
   end interface
   integer, intent(in) :: yac_instance_id
   character(len=*), intent(in) :: comp_name
   logical :: has_metadata
   has_metadata =                              &
     c_associated(                             &
       yac_cget_component_metadata_instance_c( &
         yac_instance_id, TRIM(comp_name) // c_null_char))
 end function yac_fcomponent_has_metadata_instance

 function yac_fgrid_has_metadata(grid_name) result( has_metadata )
   use, intrinsic :: iso_c_binding, only: c_null_char, c_associated
   use mo_yac_finterface, dummy => yac_fgrid_has_metadata
   use mo_yac_iso_c_helpers
   implicit none

   interface
      function yac_cget_grid_metadata_c(grid_name) result(metadata) &
           bind( c, name="yac_cget_grid_metadata")
        use, intrinsic :: iso_c_binding, only: c_char, c_ptr
        character(kind=c_char), dimension(*), intent(in) :: grid_name
        type(c_ptr) :: metadata
      end function yac_cget_grid_metadata_c
   end interface
   character(len=*), intent(in) :: grid_name
   logical :: has_metadata
   has_metadata = &
    c_associated(yac_cget_grid_metadata_c(TRIM(grid_name) // c_null_char))
 end function yac_fgrid_has_metadata

 function yac_fgrid_has_metadata_instance(yac_instance_id, grid_name) &
      result( has_metadata )
   use, intrinsic :: iso_c_binding, only: c_null_char, c_associated
   use mo_yac_finterface, dummy => yac_fgrid_has_metadata_instance
   use mo_yac_iso_c_helpers
   implicit none

   interface
      function yac_cget_grid_metadata_instance_c(yac_instance_id, grid_name) &
           result(metadata) &
           bind( c, name="yac_cget_grid_metadata_instance")
        use, intrinsic :: iso_c_binding, only: c_char, c_ptr, c_int
        integer(kind=c_int), intent(in), value :: yac_instance_id
        character(kind=c_char), dimension(*), intent(in) :: grid_name
        type(c_ptr) :: metadata
      end function yac_cget_grid_metadata_instance_c
   end interface
   integer, intent(in) :: yac_instance_id
   character(len=*), intent(in) :: grid_name
   logical :: has_metadata
   has_metadata =                         &
     c_associated(                        &
       yac_cget_grid_metadata_instance_c( &
         yac_instance_id, TRIM(grid_name) // c_null_char))
 end function yac_fgrid_has_metadata_instance

 function yac_ffield_has_metadata(comp_name, grid_name, field_name) &
      result( has_metadata )
   use, intrinsic :: iso_c_binding, only: c_null_char, c_associated
   use mo_yac_finterface, dummy => yac_ffield_has_metadata
   use mo_yac_iso_c_helpers
   implicit none

   interface
      function yac_cget_field_metadata_c(comp_name, grid_name, field_name) &
           result(metadata) &
           bind( c, name="yac_cget_field_metadata")
        use, intrinsic :: iso_c_binding, only: c_char, c_ptr
        character(kind=c_char), dimension(*), intent(in) :: comp_name
        character(kind=c_char), dimension(*), intent(in) :: grid_name
        character(kind=c_char), dimension(*), intent(in) :: field_name
        type(c_ptr) :: metadata
      end function yac_cget_field_metadata_c
   end interface
   character(len=*), intent(in) :: comp_name
   character(len=*), intent(in) :: grid_name
   character(len=*), intent(in) :: field_name
   logical :: has_metadata
   has_metadata =                        &
     c_associated(                       &
       yac_cget_field_metadata_c(        &
         TRIM(comp_name) // c_null_char, &
         TRIM(grid_name) // c_null_char, &
         TRIM(field_name) // c_null_char))
 end function yac_ffield_has_metadata

 function yac_ffield_has_metadata_instance(              &
      yac_instance_id, comp_name, grid_name, field_name) &
      result( has_metadata )
   use, intrinsic :: iso_c_binding, only: c_null_char, c_associated
   use mo_yac_finterface, dummy => yac_ffield_has_metadata_instance
   use mo_yac_iso_c_helpers
   implicit none

   interface
      function yac_cget_field_metadata_instance_c(yac_instance_id, &
                                                  comp_name,       &
                                                  grid_name,       &
                                                  field_name)      &
           result(metadata)                                        &
           bind( c, name="yac_cget_field_metadata_instance")
        use, intrinsic :: iso_c_binding, only: c_char, c_ptr, c_int
        integer(kind=c_int), intent(in), value :: yac_instance_id
        character(kind=c_char), dimension(*), intent(in) :: comp_name
        character(kind=c_char), dimension(*), intent(in) :: grid_name
        character(kind=c_char), dimension(*), intent(in) :: field_name
        type(c_ptr) :: metadata
      end function yac_cget_field_metadata_instance_c
   end interface
   integer, intent(in) :: yac_instance_id
   character(len=*), intent(in) :: comp_name
   character(len=*), intent(in) :: grid_name
   character(len=*), intent(in) :: field_name
   logical :: has_metadata
   has_metadata =                          &
     c_associated(                         &
       yac_cget_field_metadata_instance_c( &
         yac_instance_id,                  &
         TRIM(comp_name) // c_null_char,   &
         TRIM(grid_name) // c_null_char,   &
         TRIM(field_name) // c_null_char))
 end function yac_ffield_has_metadata_instance

 function yac_fget_component_metadata(comp_name) result( metadata )
   use, intrinsic :: iso_c_binding, only: c_null_char, c_ptr, c_associated
   use mo_yac_finterface, dummy => yac_fget_component_metadata
   use mo_yac_iso_c_helpers
   implicit none

   interface
      function yac_cget_component_metadata_c(comp_name) result(metadata) &
           bind( c, name="yac_cget_component_metadata")
        use, intrinsic :: iso_c_binding, only: c_char, c_ptr
        character(kind=c_char), dimension(*), intent(in) :: comp_name
        type(c_ptr) :: metadata
      end function yac_cget_component_metadata_c
   end interface
   character(len=*), intent(in) :: comp_name
   character(len=:), allocatable :: metadata
   type(c_ptr) :: c_metadata
   c_metadata = yac_cget_component_metadata_c(TRIM(comp_name) // c_null_char)
   YAC_FASSERT(c_associated(c_metadata), "ERROR(yac_fget_component_metadata): no metadata defined for component " // TRIM(comp_name))
   metadata = yac_internal_cptr2char(c_metadata)
 end function yac_fget_component_metadata

 function yac_fget_component_metadata_instance(yac_instance_id, comp_name) &
      result( metadata )

   use, intrinsic :: iso_c_binding, only: c_null_char, c_ptr, c_associated
   use mo_yac_finterface, dummy => yac_fget_component_metadata_instance
   use mo_yac_iso_c_helpers
   implicit none

   interface
      function yac_cget_component_metadata_instance_c(yac_instance_id, &
                                                      comp_name)       &
          result(metadata)                                             &
           bind( c, name="yac_cget_component_metadata_instance")
        use, intrinsic :: iso_c_binding, only: c_char, c_ptr, c_int
        integer(kind=c_int), intent(in), value :: yac_instance_id
        character(kind=c_char), dimension(*), intent(in) :: comp_name
        type(c_ptr) :: metadata
      end function yac_cget_component_metadata_instance_c
   end interface
   integer, intent(in) :: yac_instance_id
   character(len=*), intent(in) :: comp_name
   character(len=:), allocatable :: metadata
   type(c_ptr) :: c_metadata
   c_metadata =                              &
     yac_cget_component_metadata_instance_c( &
       yac_instance_id, TRIM(comp_name) // c_null_char)
   YAC_FASSERT(c_associated(c_metadata), "ERROR(yac_fget_component_metadata_instance): no metadata defined for component " // TRIM(comp_name))
   metadata = yac_internal_cptr2char(c_metadata)
 end function yac_fget_component_metadata_instance

 function yac_fget_grid_metadata(grid_name) result( metadata )
   use, intrinsic :: iso_c_binding, only: c_null_char, c_ptr, c_associated
   use mo_yac_finterface, dummy => yac_fget_grid_metadata
   use mo_yac_iso_c_helpers
   implicit none

   interface
      function yac_cget_grid_metadata_c(grid_name) result(metadata) &
           bind( c, name="yac_cget_grid_metadata")
        use, intrinsic :: iso_c_binding, only: c_char, c_ptr
        character(kind=c_char), dimension(*), intent(in) :: grid_name
        type(c_ptr) :: metadata
      end function yac_cget_grid_metadata_c
   end interface
   character(len=*), intent(in) :: grid_name
   character(len=:), allocatable :: metadata
   type(c_ptr) :: c_metadata
   c_metadata = yac_cget_grid_metadata_c(TRIM(grid_name) // c_null_char)
   YAC_FASSERT(c_associated(c_metadata), "ERROR(yac_fget_grid_metadata): no metadata defined for grid " // TRIM(grid_name))
   metadata = yac_internal_cptr2char(c_metadata)
 end function yac_fget_grid_metadata

 function yac_fget_grid_metadata_instance(yac_instance_id, grid_name) &
      result( metadata )

   use, intrinsic :: iso_c_binding, only: c_null_char, c_ptr, c_associated
   use mo_yac_finterface, dummy => yac_fget_grid_metadata_instance
   use mo_yac_iso_c_helpers
   implicit none

   interface
      function yac_cget_grid_metadata_instance_c(yac_instance_id, &
                                                 grid_name)       &
           result(metadata)                                       &
           bind( c, name="yac_cget_grid_metadata_instance")
        use, intrinsic :: iso_c_binding, only: c_char, c_ptr, c_int
        integer(kind=c_int), intent(in), value :: yac_instance_id
        character(kind=c_char), dimension(*), intent(in) :: grid_name
        type(c_ptr) :: metadata
      end function yac_cget_grid_metadata_instance_c
   end interface
   integer, intent(in) :: yac_instance_id
   character(len=*), intent(in) :: grid_name
   character(len=:), allocatable :: metadata
   type(c_ptr) :: c_metadata
   c_metadata =                         &
     yac_cget_grid_metadata_instance_c( &
       yac_instance_id, TRIM(grid_name) // c_null_char)
   YAC_FASSERT(c_associated(c_metadata), "ERROR(yac_fget_grid_metadata_instance): no metadata defined for grid " // TRIM(grid_name))
   metadata = yac_internal_cptr2char(c_metadata)
 end function yac_fget_grid_metadata_instance

 function yac_fget_field_metadata(comp_name, grid_name, field_name) &
      result( metadata )
   use, intrinsic :: iso_c_binding, only: c_null_char, c_ptr, c_associated
   use mo_yac_finterface, dummy => yac_fget_field_metadata
   use mo_yac_iso_c_helpers
   implicit none

   interface
      function yac_cget_field_metadata_c(comp_name, grid_name, field_name) &
           result(metadata) &
           bind( c, name="yac_cget_field_metadata")
        use, intrinsic :: iso_c_binding, only: c_char, c_ptr
        character(kind=c_char), dimension(*), intent(in) :: comp_name
        character(kind=c_char), dimension(*), intent(in) :: grid_name
        character(kind=c_char), dimension(*), intent(in) :: field_name
        type(c_ptr) :: metadata
      end function yac_cget_field_metadata_c
   end interface
   character(len=*), intent(in) :: comp_name
   character(len=*), intent(in) :: grid_name
   character(len=*), intent(in) :: field_name
   character(len=:), allocatable :: metadata
   type(c_ptr) :: c_metadata
   c_metadata =                         &
    yac_cget_field_metadata_c(          &
        TRIM(comp_name) // c_null_char, &
        TRIM(grid_name) // c_null_char, &
        TRIM(field_name) // c_null_char)
   YAC_FASSERT(c_associated(c_metadata), "ERROR(yac_fget_field_metadata): no metadata defined for field " // TRIM(comp_name) // "::" // TRIM(grid_name) // "::" // TRIM(field_name))
   metadata = yac_internal_cptr2char(c_metadata)
 end function yac_fget_field_metadata

 function yac_fget_field_metadata_instance(yac_instance_id, comp_name, &
      grid_name, field_name) result( metadata )

   use, intrinsic :: iso_c_binding, only: c_null_char, c_ptr, c_associated
   use mo_yac_finterface, dummy => yac_fget_field_metadata_instance
   use mo_yac_iso_c_helpers
   implicit none

   interface
      function yac_cget_field_metadata_instance_c(yac_instance_id, &
                                                  comp_name,       &
                                                  grid_name,       &
                                                  field_name)      &
           result(metadata)                                        &
           bind( c, name="yac_cget_field_metadata_instance")
        use, intrinsic :: iso_c_binding, only: c_char, c_ptr, c_int
        integer(kind=c_int), intent(in), value :: yac_instance_id
        character(kind=c_char), dimension(*), intent(in) :: comp_name
        character(kind=c_char), dimension(*), intent(in) :: grid_name
        character(kind=c_char), dimension(*), intent(in) :: field_name
        type(c_ptr) :: metadata
      end function yac_cget_field_metadata_instance_c
   end interface
   integer, intent(in) :: yac_instance_id
   character(len=*), intent(in) :: comp_name
   character(len=*), intent(in) :: grid_name
   character(len=*), intent(in) :: field_name
   character(len=:), allocatable :: metadata
   type(c_ptr) :: c_metadata
   c_metadata =                         &
    yac_cget_field_metadata_instance_c( &
        yac_instance_id,                &
        TRIM(comp_name) // c_null_char, &
        TRIM(grid_name) // c_null_char, &
        TRIM(field_name) // c_null_char)
   YAC_FASSERT(c_associated(c_metadata), "ERROR(yac_fget_field_metadata_instance): no metadata defined for field " // TRIM(comp_name) // "::" // TRIM(grid_name) // "::" // TRIM(field_name))
   metadata = yac_internal_cptr2char(c_metadata)
 end function yac_fget_field_metadata_instance

 ! ---------------------------------------------------------------------

function yac_fget_start_datetime () result (start_datetime_string)

  use, intrinsic :: iso_c_binding, only : c_ptr
  use mo_yac_iso_c_helpers
  use mo_yac_finterface, dummy => yac_fget_start_datetime

  implicit none

  interface
     function yac_cget_start_datetime_c() &
        bind ( c, name='yac_cget_start_datetime' )

     use, intrinsic :: iso_c_binding, only : c_ptr
     type(c_ptr) :: yac_cget_start_datetime_c

     end function yac_cget_start_datetime_c
  end interface

  type (c_ptr)                   :: c_string_ptr
  character (len=:), ALLOCATABLE :: start_datetime_string

  c_string_ptr = yac_cget_start_datetime_c()
  start_datetime_string = yac_internal_cptr2char(c_string_ptr)

end function yac_fget_start_datetime

function yac_fget_start_datetime_instance (yac_instance_id) &
  result (start_datetime_string)

  use, intrinsic :: iso_c_binding, only : c_ptr
  use mo_yac_iso_c_helpers
  use mo_yac_finterface, dummy => yac_fget_start_datetime_instance

  implicit none

  interface
     function yac_cget_start_datetime_instance_c(yac_instance_id) &
        bind ( c, name='yac_cget_start_datetime_instance' )

     use, intrinsic :: iso_c_binding, only : c_ptr, c_int
     integer ( kind=c_int ), value :: yac_instance_id
     type(c_ptr)                   :: yac_cget_start_datetime_instance_c

     end function yac_cget_start_datetime_instance_c
  end interface

  integer, intent(in) :: yac_instance_id !< [IN]  YAC instance identifier

  type (c_ptr)                   :: c_string_ptr
  character (len=:), ALLOCATABLE :: start_datetime_string

  c_string_ptr = yac_cget_start_datetime_instance_c(yac_instance_id)
  start_datetime_string = yac_internal_cptr2char(c_string_ptr)

end function yac_fget_start_datetime_instance

! ---------------------------------------------------------------------

function yac_fget_end_datetime () result (end_datetime_string)

  use, intrinsic :: iso_c_binding, only : c_ptr
  use mo_yac_iso_c_helpers
  use mo_yac_finterface, dummy => yac_fget_end_datetime

  implicit none

  interface
     function yac_cget_end_datetime_c() &
        bind ( c, name='yac_cget_end_datetime' )

     use, intrinsic :: iso_c_binding, only : c_ptr
     type(c_ptr) :: yac_cget_end_datetime_c

     end function yac_cget_end_datetime_c
  end interface

  type (c_ptr)                   :: c_string_ptr
  character (len=:), ALLOCATABLE :: end_datetime_string

  c_string_ptr = yac_cget_end_datetime_c()
  end_datetime_string = yac_internal_cptr2char(c_string_ptr)

end function yac_fget_end_datetime

function yac_fget_end_datetime_instance (yac_instance_id) &
  result (end_datetime_string)

  use, intrinsic :: iso_c_binding, only : c_ptr
  use mo_yac_iso_c_helpers
  use mo_yac_finterface, dummy => yac_fget_end_datetime_instance

  implicit none

  interface
     function yac_cget_end_datetime_instance_c(yac_instance_id) &
        bind ( c, name='yac_cget_end_datetime_instance' )

     use, intrinsic :: iso_c_binding, only : c_ptr, c_int
     integer ( kind=c_int ), value :: yac_instance_id
     type(c_ptr)                   :: yac_cget_end_datetime_instance_c

     end function yac_cget_end_datetime_instance_c
  end interface

  integer, intent(in) :: yac_instance_id !< [IN]  YAC instance identifier

  type (c_ptr)                   :: c_string_ptr
  character (len=:), ALLOCATABLE :: end_datetime_string

  c_string_ptr = yac_cget_end_datetime_instance_c(yac_instance_id)
  end_datetime_string = yac_internal_cptr2char(c_string_ptr)

end function yac_fget_end_datetime_instance

! ---------------------------------------------------------------------

subroutine yac_redirstdout_f ( filestem, parallel, my_pe, npes )

  use mo_yac_finterface, dummy => yac_redirstdout_f
  use, intrinsic :: iso_c_binding, only : c_null_char
  use mo_yac_iso_c_helpers

  implicit none

  interface
     subroutine yac_redirstdout_c ( filestem, parallel, my_pe, npes ) &
    &      bind ( c, name='yac_redirstdout' )

       use, intrinsic :: iso_c_binding, only : c_int, c_char

       character ( kind=c_char ), dimension(*) :: filestem
       integer ( kind=c_int ), value           :: parallel
       integer ( kind=c_int ), value           :: my_pe
       integer ( kind=c_int ), value           :: npes

     end subroutine yac_redirstdout_c
  end interface

  character(len=*), intent(in) :: filestem
  integer, intent(in) :: parallel
  integer, intent(in) :: my_pe
  integer, intent(in) :: npes

  call yac_check_strlength ( filestem )

  call yac_redirstdout_c ( &
     TRIM(filestem) // c_null_char, parallel, my_pe, npes)

end subroutine yac_redirstdout_f

! ---------------------------------------------------------------------

subroutine yac_fget_interp_stack_config ( interp_stack_config_id )

  use mo_yac_finterface, dummy => yac_fget_interp_stack_config

  implicit none

  interface

      subroutine yac_cget_interp_stack_config_c ( interp_stack_config_id ) &
        bind ( c, name='yac_cget_interp_stack_config' )
        use, intrinsic :: iso_c_binding, only : c_int

        integer ( kind=c_int ) :: interp_stack_config_id

      end subroutine yac_cget_interp_stack_config_c

   end interface

  integer, intent(out) :: interp_stack_config_id

  call yac_cget_interp_stack_config_c ( interp_stack_config_id )

end subroutine yac_fget_interp_stack_config

subroutine yac_ffree_interp_stack_config ( interp_stack_config_id )

  use mo_yac_finterface, dummy => yac_ffree_interp_stack_config

  implicit none

  interface

      subroutine yac_cfree_interp_stack_config_c ( interp_stack_config_id ) &
        bind ( c, name='yac_cfree_interp_stack_config' )
        use, intrinsic :: iso_c_binding, only : c_int

        integer ( kind=c_int ), value :: interp_stack_config_id

      end subroutine yac_cfree_interp_stack_config_c

   end interface

  integer, intent(in) :: interp_stack_config_id

  call yac_cfree_interp_stack_config_c ( interp_stack_config_id )

end subroutine yac_ffree_interp_stack_config

subroutine yac_fadd_interp_stack_config_average ( interp_stack_config_id, &
     reduction_type, partial_coverage)

  use mo_yac_finterface, dummy => yac_fadd_interp_stack_config_average

  implicit none

  interface

     subroutine yac_cadd_interp_stack_config_average_c (             &
          interp_stack_config_id, reduction_type, partial_coverage ) &
        bind ( c, name='yac_cadd_interp_stack_config_average' )
        use, intrinsic :: iso_c_binding, only : c_int

        integer ( kind=c_int ), value :: interp_stack_config_id
        integer ( kind=c_int ), value :: reduction_type
        integer ( kind=c_int ), value :: partial_coverage

      end subroutine yac_cadd_interp_stack_config_average_c

   end interface

   integer, intent(in)          :: interp_stack_config_id
   integer, intent(in)          :: reduction_type
   integer, intent(in)          :: partial_coverage

  call yac_cadd_interp_stack_config_average_c ( &
         interp_stack_config_id, reduction_type, partial_coverage )

end subroutine yac_fadd_interp_stack_config_average

subroutine yac_fadd_interp_stack_config_nnn(interp_stack_config_id, &
     type, n, scale)

  use mo_yac_finterface, dummy => yac_fadd_interp_stack_config_nnn

  implicit none

  interface

     subroutine yac_cadd_interp_stack_config_nnn_c( &
          interp_stack_config_id, type, n, scale)   &
        bind ( c, name='yac_cadd_interp_stack_config_nnn' )
        use, intrinsic :: iso_c_binding, only : c_int, c_double

        integer ( kind=c_int ), value :: interp_stack_config_id
        integer ( kind=c_int ), value :: type
        integer ( kind=c_int ), value :: n
        real ( kind=c_double ), value :: scale

      end subroutine yac_cadd_interp_stack_config_nnn_c

   end interface

   integer, intent(in)          :: interp_stack_config_id
   integer, intent(in)          :: type
   integer, intent(in)          :: n
   double precision, intent(in) :: scale

  call yac_cadd_interp_stack_config_nnn_c( &
         interp_stack_config_id, type, n, scale)

end subroutine yac_fadd_interp_stack_config_nnn

subroutine yac_fadd_interp_stack_config_conservative( &
     interp_stack_config_id, order, enforced_conserv, &
     partial_coverage, normalization)

  use mo_yac_finterface, dummy => yac_fadd_interp_stack_config_conservative

  implicit none

  interface

     subroutine yac_cadd_interp_stack_config_conservative_c( &
         interp_stack_config_id, order, enforced_conserv,    &
         partial_coverage, normalization)                    &
        bind ( c, name='yac_cadd_interp_stack_config_conservative' )
        use, intrinsic :: iso_c_binding, only : c_int, c_double

        integer ( kind=c_int ), value :: interp_stack_config_id
        integer ( kind=c_int ), value :: order
        integer ( kind=c_int ), value :: enforced_conserv
        integer ( kind=c_int ), value :: partial_coverage
        integer ( kind=c_int ), value :: normalization

      end subroutine yac_cadd_interp_stack_config_conservative_c

   end interface

   integer, intent(in)     :: interp_stack_config_id
   integer, intent(in)     :: order
   integer, intent(in)     :: enforced_conserv
   integer, intent(in)     :: partial_coverage
   integer, intent(in)     :: normalization

  call yac_cadd_interp_stack_config_conservative_c(       &
         interp_stack_config_id, order, enforced_conserv, &
         partial_coverage, normalization)

end subroutine yac_fadd_interp_stack_config_conservative

subroutine yac_fadd_interp_stack_config_spmap(interp_stack_config_id, &
     spread_distance, max_search_distance, weight_type)

  use mo_yac_finterface, dummy => yac_fadd_interp_stack_config_spmap

  implicit none

  interface

     subroutine yac_cadd_interp_stack_config_spmap_c( &
          interp_stack_config_id, spread_distance,    &
          max_search_distance, weight_type)           &
          bind ( c, name='yac_cadd_interp_stack_config_spmap' )
       use, intrinsic :: iso_c_binding, only : c_int, c_double

       integer ( kind=c_int ), value :: interp_stack_config_id
       real ( kind=c_double ), value :: spread_distance
       real ( kind=c_double ), value :: max_search_distance
       integer ( kind=c_int ), value :: weight_type

     end subroutine yac_cadd_interp_stack_config_spmap_c

  end interface

  integer, intent(in)          :: interp_stack_config_id
  double precision, intent(in) :: spread_distance
  double precision, intent(in) :: max_search_distance
  integer, intent(in)          :: weight_type

  call yac_cadd_interp_stack_config_spmap_c(      &
         interp_stack_config_id, spread_distance, &
         max_search_distance, weight_type)

end subroutine yac_fadd_interp_stack_config_spmap

subroutine yac_fadd_interp_stack_config_hcsbb(interp_stack_config_id)

  use mo_yac_finterface, dummy => yac_fadd_interp_stack_config_hcsbb

  implicit none

  interface

     subroutine yac_cadd_interp_stack_config_hcsbb_c( &
          interp_stack_config_id)                     &
          bind ( c, name='yac_cadd_interp_stack_config_hcsbb' )
       use, intrinsic :: iso_c_binding, only : c_int

       integer ( kind=c_int ), value :: interp_stack_config_id

     end subroutine yac_cadd_interp_stack_config_hcsbb_c

  end interface

  integer, intent(in) :: interp_stack_config_id

  call yac_cadd_interp_stack_config_hcsbb_c(interp_stack_config_id)

end subroutine yac_fadd_interp_stack_config_hcsbb

subroutine yac_fadd_interp_stack_config_user_file(interp_stack_config_id, &
     filename, src_grid_name, tgt_grid_name)

  use mo_yac_finterface, dummy => yac_fadd_interp_stack_config_user_file
  use, intrinsic :: iso_c_binding, only : c_null_char

  implicit none

  interface

     subroutine yac_cadd_interp_stack_config_user_file_c(  &
          interp_stack_config_id, filename, src_grid_name, &
          tgt_grid_name)                                   &
          bind ( c, name='yac_cadd_interp_stack_config_user_file' )
       use, intrinsic :: iso_c_binding, only : c_int, c_char

       integer ( kind=c_int ), value :: interp_stack_config_id
       character (kind=c_char), dimension(*) :: filename
       character (kind=c_char), dimension(*) :: src_grid_name
       character (kind=c_char), dimension(*) :: tgt_grid_name

     end subroutine yac_cadd_interp_stack_config_user_file_c

  end interface

  integer, intent(in)           :: interp_stack_config_id
  character (len=*), intent(in) :: filename
  character (len=*), intent(in) :: src_grid_name
  character (len=*), intent(in) :: tgt_grid_name

  call yac_cadd_interp_stack_config_user_file_c( &
         interp_stack_config_id,                 &
         TRIM(filename) // c_null_char,          &
         TRIM(src_grid_name) // c_null_char,     &
         TRIM(tgt_grid_name) // c_null_char)

end subroutine yac_fadd_interp_stack_config_user_file

subroutine yac_fadd_interp_stack_config_fixed(interp_stack_config_id, &
     val)

  use mo_yac_finterface, dummy => yac_fadd_interp_stack_config_fixed

  implicit none

  interface

     subroutine yac_cadd_interp_stack_config_fixed_c( &
          interp_stack_config_id, val)                &
          bind ( c, name='yac_cadd_interp_stack_config_fixed' )
       use, intrinsic :: iso_c_binding, only : c_int, c_double

       integer ( kind=c_int ), value :: interp_stack_config_id
       real ( kind=c_double ), value :: val

     end subroutine yac_cadd_interp_stack_config_fixed_c

  end interface

  integer, intent(in)          :: interp_stack_config_id
  double precision, intent(in) :: val

  call yac_cadd_interp_stack_config_fixed_c( &
         interp_stack_config_id, val)

end subroutine yac_fadd_interp_stack_config_fixed

subroutine yac_fadd_interp_stack_config_creep(interp_stack_config_id, &
     creep_distance)

  use mo_yac_finterface, dummy => yac_fadd_interp_stack_config_creep

  implicit none

  interface

     subroutine yac_cadd_interp_stack_config_creep_c( &
          interp_stack_config_id, creep_distance)     &
          bind ( c, name='yac_cadd_interp_stack_config_creep' )
       use, intrinsic :: iso_c_binding, only : c_int

       integer ( kind=c_int ), value :: interp_stack_config_id
       integer ( kind=c_int ), value :: creep_distance

     end subroutine yac_cadd_interp_stack_config_creep_c

  end interface

  integer, intent(in) :: interp_stack_config_id
  integer, intent(in) :: creep_distance

  call yac_cadd_interp_stack_config_creep_c( &
         interp_stack_config_id, creep_distance)

end subroutine yac_fadd_interp_stack_config_creep
