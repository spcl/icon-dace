! mo_machine.f90 - Machine related information
!
! Copyright (C) 2014, MPI-M
! SPDX-License-Identifier: BSD-3-Clause
! See ./LICENSES/ for license information
!_________________________________________

#ifdef __xlC__
@PROCESS STRICT
#endif
MODULE mo_machine
  !
  !  This routine originates (year 2014) from MPI-ESM, the Earth System Model of the 
  !  Max Planck Institute for Meteorology (Mauritsen et al. 2019). 
  !  Reference: Mauritsen, T., et al. (2019) Developments in the MPI-M Earth System Model 
  !  version 1.2 (MPI-ESM1.2) and its response to increasing CO2. J. Adv. Model. Earth Syst., 11, 
  !  doi: 10.1029/2018MS001400.

  USE mo_kind,      ONLY: sp, dp, i4, i8
  USE mo_mpi,       ONLY: p_parallel_io
  USE mo_exception, ONLY: message_text, message, finish

  IMPLICIT NONE

  PUBLIC

  INTEGER :: mp_integer, mp_int4, mp_int8
  INTEGER :: mp_real, mp_real4, mp_real8
  INTEGER :: mp_logical

  REAL (dp) :: epsi     ! relative machine precision
  REAL (dp) :: sfmin    ! safe minimum, such that 1/sfmin does not overflow
  INTEGER   :: base     ! base of the machine
  REAL (dp) :: prec     ! epsi*base
  INTEGER   :: mantissa ! number of (base) digits in the mantissa
  INTEGER   :: emin     ! minimum exponent before (gradual) underflow
  REAL (dp) :: rmin     ! underflow threshold - base**(emin-1)
  INTEGER   :: emax     ! largest exponent before overflow
  REAL (dp) :: rmax     ! overflow threshold  - (base**emax)*(1-epsi)

  LOGICAL, PRIVATE :: lrnd, lieee 

CONTAINS

  SUBROUTINE machine_setup

    INTEGER      :: io_size
    INTEGER      :: integer_byte_size, integer_io_size
    INTEGER (i4) :: int4_bitsize
    INTEGER (i8) :: int8_bitsize

    INTEGER      :: ii   = 0      
    INTEGER (i4) :: int4 = 0_i4
    INTEGER (i8) :: int8 = 0_i8

    REAL         :: rr = 0.0
    REAL (sp)    :: ss = 0.0_sp
    REAL (dp)    :: dd = 0.0_dp

    LOGICAL      :: ll = .TRUE.

    integer_byte_size = BIT_SIZE(ii)/8

    INQUIRE (iolength=io_size) ii
    integer_io_size = io_size
    mp_integer = io_size/integer_io_size*integer_byte_size

    INQUIRE (iolength=io_size) int4
    mp_int4 = io_size/integer_io_size*integer_byte_size
    INQUIRE (iolength=io_size) int8
    mp_int8 = io_size/integer_io_size*integer_byte_size

    INQUIRE (iolength=io_size) rr
    mp_real = io_size/integer_io_size*integer_byte_size

    INQUIRE (iolength=io_size) ss
    mp_real4 = io_size/integer_io_size*integer_byte_size
    INQUIRE (iolength=io_size) dd
    mp_real8 = io_size/integer_io_size*integer_byte_size

    INQUIRE (iolength=io_size) ll
    mp_logical = io_size/integer_io_size*integer_byte_size

    int4_bitsize = BIT_SIZE(int4)
    int8_bitsize = BIT_SIZE(int8)

    IF (int4_bitsize /= 8*mp_int4) THEN
       IF (p_parallel_io) THEN
          CALL message('',' There is a problem with the compiler options/setup or it is an very old CRAY!')
          WRITE (message_text,'(a,i0,a,i0,a)') &
               ' The bit size of an INTEGER (4 byte range) is ', &
               int4_bitsize, ' iolength returns ', 8*mp_int4, 'bit.'
          CALL message('',TRIM(message_text))
       END IF
    END IF

    IF (int8_bitsize /= 8*mp_int8) THEN
       IF (p_parallel_io) THEN
          CALL message('',' There is a serious problem with the compiler options/setup!')
          WRITE (message_text,'(a,i0,a,i0,a)') &
               ' The bit size of an INTEGER (8 byte range) is ', &  
               int8_bitsize, ' iolength returns ', 8*mp_int8, 'bit.'
          CALL message('',TRIM(message_text))
       END IF
    END IF

    IF (p_parallel_io) THEN
       CALL message('',' Automatic determined byte size of types for I/O:')
       WRITE (message_text, '(a,i3)') '  INTEGER                : ', mp_integer
       CALL message('',TRIM(message_text))
       WRITE (message_text, '(a,i3)') '  INTEGER (4 byte range) : ', mp_int4
       CALL message('',TRIM(message_text))
       WRITE (message_text, '(a,i3)') '  INTEGER (8 byte range) : ', mp_int8
       CALL message('',TRIM(message_text))
       WRITE (message_text, '(a,i3)') '  REAL                   : ', mp_real  
       CALL message('',TRIM(message_text))
       WRITE (message_text, '(a,i3)') '  REAL (single)          : ', mp_real4  
       CALL message('',TRIM(message_text))
       WRITE (message_text, '(a,i3)') '  REAL (double)          : ', mp_real8
       CALL message('',TRIM(message_text))
       WRITE (message_text, '(a,i3)') '  LOGICAL                : ', mp_logical
       CALL message('',TRIM(message_text))
       CALL message('','')
    END IF

    IF (mp_real8 < 8) THEN
       IF (p_parallel_io) THEN
          CALL message('',' There is a serious problem with the compiler options/ setup!')
          WRITE (message_text,'(a,i0)') ' The byte size of an REAL(kind=dp) is only ', mp_real8
          CALL message('',TRIM(message_text))
       END IF
       CALL finish('machine_setup','setup error.')
    END IF

    sfmin    = SPACING (dd)
    base     = RADIX (dd)
    prec     = EPSILON (dd)
    mantissa = DIGITS (dd)
    emin     = MINEXPONENT (dd)
    rmin     = TINY (dd)
    emax     = MAXEXPONENT (dd)
    rmax     = HUGE (dd)

    CALL rounding

    IF (lrnd) THEN
       epsi = (REAL(base,dp)**(1-mantissa))/2
    ELSE
       epsi = REAL(base,dp)**(1-mantissa)
    END IF

    IF (p_parallel_io) THEN
       CALL message('',' Automatic determined machine constants (8 byte REAL):')
       WRITE (message_text,'(a,e25.15)')  '  Epsilon                      = ', epsi
       CALL message('',TRIM(message_text))
       WRITE (message_text,'(a,e25.15)')  '  Precision                    = ', prec
       CALL message('',TRIM(message_text))
       WRITE (message_text,'(a,e25.15)')  '  Safe minimum                 = ', sfmin
       CALL message('',TRIM(message_text))
       WRITE (message_text,'(a,e25.15)')  '  Reciprocal of safe minimum   = ', 1/sfmin
       CALL message('',TRIM(message_text))
       WRITE (message_text,'(a,e25.15)')  '  Underflow threshold          = ', rmin
       CALL message('',TRIM(message_text))
       WRITE (message_text,'(a,e25.15)')  '  Overflow threshold           = ', rmax
       CALL message('','')
       WRITE (message_text,'(a,i6)')      '  Number of digits in mantissa = ', mantissa
       CALL message('',TRIM(message_text))
       WRITE (message_text,'(a,i6)')      '  Base                         = ', base
       CALL message('',TRIM(message_text))
       WRITE (message_text,'(a,i6)')      '  Minimum exponent             = ', emin
       CALL message('',TRIM(message_text))
       WRITE (message_text,'(a,i6)')      '  Largest exponent             = ', emax
       CALL message('',TRIM(message_text))

       CALL message('','')

       IF (lieee) THEN
          CALL message('',' Rounding mode according the IEEE definition.')
       ELSE
          IF (lrnd) THEN
             CALL message('','Rounding mode gives proper rounding.')
          ELSE
             CALL message('','Rounding mode is chopping.')
          END IF
       END IF
       CALL message('','')
    END IF

  END SUBROUTINE machine_setup

  ! add_numbers necessary to prevent wrong results due to optimization

  FUNCTION add_numbers (a, b) RESULT (c)

    REAL (dp), INTENT(IN) :: a, b
    REAL (dp)             :: c  

    c = a+b

  END FUNCTION add_numbers

  SUBROUTINE rounding 

    ! lrnd    (output) LOGICAL
    !         Specifies whether proper rounding  ( lrnd = .TRUE. )  or
    !         chopping  ( lrnd = .FALSE. )  occurs in addition. This may not
    !         be a reliable guide to the way in which the machine performs
    !         its arithmetic.
    !
    ! lieee   (output) LOGICAL
    !         Specifies whether rounding appears to be done in the IEEE
    !         'round to nearest' style.

    REAL (dp) :: one = 1.0_dp
    REAL (dp) :: a   = 1.0_dp
    REAL (dp) :: b   = 1.0_dp
    REAL (dp) :: c   = 1.0_dp

    REAL (dp) :: qtr, savec, f, t1, t2
    INTEGER  :: ibeta


    ! Throughout this routine  we use the function add_numbers to ensure
    ! that relevant values are  stored and not held in registers,  or
    ! are not affected by optimizers.

    ! Compute  a = 2.0**m  with the  smallest positive integer m such
    ! that
    ! 
    ! fl( a + 1.0_dp ) = a.
    ! 

    DO WHILE (c == 1)
       a = 2*a
       c = add_numbers (a, one)
       c = add_numbers (c, -a)
    END DO

    ! Now compute  b = 2.0**m  with the smallest positive integer m
    ! such that
    ! 
    ! fl( a + b ) .gt. a. 

    c = add_numbers (a, b)
    DO WHILE (c == a) 
       b = 2*b
       c = add_numbers (a, b)
    END DO

    ! Now compute the base.  a and c are neighbouring floating point
    ! numbers  in the  interval (ibeta**t, ibeta**( t + 1 ))  and so
    ! their difference is beta. Adding 0.25 to c is to ensure that it
    ! is truncated to beta and not ( ibeta - 1 ).

    qtr = one/4
    savec = c
    c = add_numbers (c, -a)
    ibeta = c+qtr

    ! Now determine whether rounding or chopping occurs,  by adding a
    ! bit  less  than  ibeta/2  and a  bit  more  than  ibeta/2  to  a.

    b = ibeta
    f = add_numbers (b/2, -b/100)
    c = add_numbers (f, a)
    IF (c == a) THEN
       lrnd = .TRUE.
    ELSE
       lrnd = .FALSE.
    END IF
    f = add_numbers (b/2, b/100)
    c = add_numbers (f, a)
    IF ((lrnd) .AND. (c == a)) lrnd = .FALSE.

    ! Try and decide whether rounding is done in the IEEE round to
    ! nearest style. b/2 is half a unit in the last place of the two
    ! numbers a and savec. Furthermore, a is even, i.e. has last  bit
    ! zero, and savec is odd. Thus adding b/2 to a should not  change
    ! a, but adding b/2 to savec should change savec.

    t1 = add_numbers (b/2, a)
    t2 = add_numbers (b/2, savec)
    lieee = (t1 == a) .AND. (t2 > savec) .AND. lrnd

  END SUBROUTINE rounding

END MODULE mo_machine
