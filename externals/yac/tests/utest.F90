!!
!! @file test_quicksort.c
!!
!! @copyright Copyright  (C)  2013 DKRZ, MPI-M
!!
!! @author Moritz Hanke <hanke@dkrz.de>
!!         Rene Redler  <rene.redler@mpimet.mpg.de>
!!
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


module utest

   implicit none

   private

   public :: start_test, stop_test, start_timer, stop_timer, &
             test_v_r, test_r, exit_tests

   integer :: utest_testcount, utest_errorcount
   real    :: start, end

contains

   subroutine start_test (name)
      character(len=*), intent(in) :: name
      utest_errorcount = 0
      utest_testcount = 0
#ifdef VERBOSE
      if (.true.) then
#else
      if (.false.) then
#endif
         print *, "testing ", name, ":"
      endif
   end subroutine start_test

   subroutine stop_test ()
#ifdef VERBOSE
      if (.true.) then
#else
      if (.false.) then
#endif
         if (utest_errorcount == 0) then
            print *, "all",utest_testcount, "test(s) succeeded"
         else
            write(0,*) utest_errorcount, " out of ", &
                       utest_testcount, " test(s) failed"
         endif
      endif
   end subroutine stop_test

   subroutine start_timer (name)
      character(len=*), intent(in) :: name
      call cpu_time(start)
#ifdef VERBOSE
      if (.true.) then
#else
      if (.false.) then
#endif
         print *, "start timing ", name
      endif
   end subroutine start_timer

   subroutine stop_timer ()
      call cpu_time(end)
#ifdef VERBOSE
      if (.true.) then
#else
      if (.false.) then
#endif
         print '(" elapsed time:",f6.3," seconds.")',end-start
      endif
   end subroutine stop_timer

   subroutine test_v_r(arg, file, line)

      logical, intent(in) :: arg
      integer, intent(in) :: line
      character(len=*), intent(in) :: file

      utest_testcount = utest_testcount + 1
      if (arg) then
#ifdef VERBOSE
         if (.true.) then
#else
         if (.false.) then
#endif
            print *, "- test", utest_testcount, "succeeded"
         endif
         return
      else
         utest_errorcount = utest_errorcount + 1
#ifdef VERBOSE
         if (.true.) then
#else
         if (.false.) then
#endif
            write(0,*) "- test", utest_testcount, "failed (",file," :", line, ")"
         endif
      endif

   end subroutine test_v_r

   subroutine test_r(arg, file, line)

      logical, intent(in) :: arg
      integer, intent(in) :: line
      character(len=*), intent(in) :: file

      utest_testcount = utest_testcount + 1
      if (arg) then
#ifdef VERBOSE
         if (.true.) then
#else
         if (.false.) then
#endif
            print *, "- test", utest_testcount, "succeeded"
         endif
         return
      else
         utest_errorcount = utest_errorcount + 1
#ifdef VERBOSE
         if (.true.) then
#else
         if (.false.) then
#endif
            write(0,*) "- test", utest_testcount, "failed (",file," :", line, ")"
         endif
      endif

   end subroutine test_r

   subroutine exit_tests()

      if (utest_errorcount > 0) then
         stop 1
      endif
   end subroutine exit_tests

end module utest
