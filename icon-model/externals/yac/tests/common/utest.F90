! Copyright (c) 2024 The YAC Authors
!
! SPDX-License-Identifier: BSD-3-Clause

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
      print *, "testing ", name, ":"
   end subroutine start_test

   subroutine stop_test ()
      if (utest_errorcount == 0) then
         print *, "all",utest_testcount, "test(s) succeeded"
      else
         write(0,*) utest_errorcount, " out of ", &
                    utest_testcount, " test(s) failed"
      endif
   end subroutine stop_test

   subroutine start_timer (name)
      character(len=*), intent(in) :: name
      call cpu_time(start)
      print *, "start timing ", name
   end subroutine start_timer

   subroutine stop_timer ()
      call cpu_time(end)
      print '(" elapsed time:",f6.3," seconds.")',end-start
   end subroutine stop_timer

   subroutine test_v_r(arg, file, line)

      logical, intent(in) :: arg
      integer, intent(in) :: line
      character(len=*), intent(in) :: file

      utest_testcount = utest_testcount + 1
      if (arg) then
         print *, "- test", utest_testcount, "succeeded"
         return
      else
         utest_errorcount = utest_errorcount + 1
         write(0,*) "- test", utest_testcount, "failed (",file," :", line, ")"
      endif

   end subroutine test_v_r

   subroutine test_r(arg, file, line)

      logical, intent(in) :: arg
      integer, intent(in) :: line
      character(len=*), intent(in) :: file

      utest_testcount = utest_testcount + 1
      if (arg) then
         print *, "- test", utest_testcount, "succeeded"
         return
      else
         utest_errorcount = utest_errorcount + 1
         write(0,*) "- test", utest_testcount, "failed (",file," :", line, ")"
      endif

   end subroutine test_r

   subroutine exit_tests()

      if (utest_errorcount > 0) then
         stop 1
      endif
   end subroutine exit_tests

end module utest
