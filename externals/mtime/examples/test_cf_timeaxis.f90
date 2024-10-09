!! Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
!!
!! SPDX-License-Identifier: BSD-3-Clause
!!
PROGRAM test_cf_timeaxis

  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: int64, wp => real64

  USE mtime

  IMPLICIT NONE

  TYPE(juliandelta), PARAMETER :: timeunit_day = juliandelta('+', 0_int64, 86400000_int64)
  TYPE(juliandelta), PARAMETER :: timeunit_hour = juliandelta('+', 0_int64, 3600000_int64)
  TYPE(juliandelta), PARAMETER :: timeunit_minute = juliandelta('+', 0_int64, 60000_int64)
  TYPE(juliandelta), PARAMETER :: timeunit_second = juliandelta('+', 0_int64, 1000_int64)

  ! epoch: 1870-01-01T00:00
  TYPE(julianday) :: epoch = julianday(2404063, 43200000)

  ! example 1:  "days since 1870-01-01T00:00"
  !
  CHARACTER(len=4) :: time1_unit = 'days'
  REAL(wp) :: time1(148) = [38730.5, 38760.5, 38790.5, 38821.0, 38851.5, &
       &      38882.0, 38912.5, 38943.5, 38974.0, 39004.5, 39035.0, 39065.5, &
       &      39096.5, 39126.0, 39155.5, 39186.0, 39216.5, 39247.0, 39277.5, &
       &      39308.5, 39339.0, 39369.5, 39400.0, 39430.5, 39461.5, 39491.0, &
       &      39520.5, 39551.0, 39581.5, 39612.0, 39642.5, 39673.5, 39704.0, &
       &      39734.5, 39765.0, 39795.5, 39826.5, 39856.0, 39885.5, 39916.0, &
       &      39946.5, 39977.0, 40007.5, 40038.5, 40069.0, 40099.5, 40130.0, &
       &      40160.5, 40191.5, 40221.5, 40251.5, 40282.0, 40312.5, 40343.0, &
       &      40373.5, 40404.5, 40435.0, 40465.5, 40496.0, 40526.5, 40557.5, &
       &      40587.0, 40616.5, 40647.0, 40677.5, 40708.0, 40738.5, 40769.5, &
       &      40800.0, 40830.5, 40861.0, 40891.5, 40922.5, 40952.0, 40981.5, &
       &      41012.0, 41042.5, 41073.0, 41103.5, 41134.5, 41165.0, 41195.5, &
       &      41226.0, 41256.5, 41287.5, 41317.0, 41346.5, 41377.0, 41407.5, &
       &      41438.0, 41468.5, 41499.5, 41530.0, 41560.5, 41591.0, 41621.5, &
       &      41652.5, 41682.5, 41712.5, 41743.0, 41773.5, 41804.0, 41834.5, &
       &      41865.5, 41896.0, 41926.5, 41957.0, 41987.5, 42018.5, 42048.0, &
       &      42077.5, 42108.0, 42138.5, 42169.0, 42199.5, 42230.5, 42261.0, &
       &      42291.5, 42322.0, 42352.5, 42383.5, 42413.0, 42442.5, 42473.0, &
       &      42503.5, 42534.0, 42564.5, 42595.5, 42626.0, 42656.5, 42687.0, &
       &      42717.5, 42748.5, 42778.0, 42807.5, 42838.0, 42868.5, 42899.0, &
       &      42929.5, 42960.5, 42991.0, 43021.5, 43052.0, 43082.5, 43113.5, &
       &      43143.5, 43173.5, 43204.0]

  ! example 2:
  CHARACTER(len=7) :: time2_unit = 'seconds'
  REAL(wp) ::  time2(125) = [0.0, 21600.0, 43200.0, 64800.0, 86400.0, &
       &       108000.0, 129600.0, 151200.0, 172800.0, 194400.0, 216000.0, &
       &       237600.0, 259200.0, 280800.0, 302400.0, 324000.0, 345600.0, &
       &       367200.0, 388800.0, 410400.0, 432000.0, 453600.0, 475200.0, &
       &       496800.0, 518400.0, 540000.0, 561600.0, 583200.0, 604800.0, &
       &       626400.0, 648000.0, 669600.0, 691200.0, 712800.0, 734400.0, &
       &       756000.0, 777600.0, 799200.0, 820800.0, 842400.0, 864000.0, &
       &       885600.0, 907200.0, 928800.0, 950400.0, 972000.0, 993600.0, &
       &      1015200.0, 1036800.0, 1058400.0, 1080000.0, 1101600.0, 1123200.0, &
       &      1144800.0, 1166400.0, 1188000.0, 1209600.0, 1231200.0, 1252800.0, &
       &      1274400.0, 1296000.0, 1317600.0, 1339200.0, 1360800.0, 1382400.0, &
       &      1404000.0, 1425600.0, 1447200.0, 1468800.0, 1490400.0, 1512000.0, &
       &      1533600.0, 1555200.0, 1576800.0, 1598400.0, 1620000.0, 1641600.0, &
       &      1663200.0, 1684800.0, 1706400.0, 1728000.0, 1749600.0, 1771200.0, &
       &      1792800.0, 1814400.0, 1836000.0, 1857600.0, 1879200.0, 1900800.0, &
       &      1922400.0, 1944000.0, 1965600.0, 1987200.0, 2008800.0, 2030400.0, &
       &      2052000.0, 2073600.0, 2095200.0, 2116800.0, 2138400.0, 2160000.0, &
       &      2181600.0, 2203200.0, 2224800.0, 2246400.0, 2268000.0, 2289600.0, &
       &      2311200.0, 2332800.0, 2354400.0, 2376000.0, 2397600.0, 2419200.0, &
       &      2440800.0, 2462400.0, 2484000.0, 2505600.0, 2527200.0, 2548800.0, &
       &      2570400.0, 2592000.0, 2613600.0, 2635200.0, 2656800.0, 2678400.0]

  PRINT *, epoch%day, epoch%ms

  iteration1: BLOCK

    TYPE(juliandelta) :: time_multiplicator
    REAL(wp) :: s, s0, days, ms
    INTEGER :: i

    SELECT CASE (time1_unit)
    CASE ('days')
      time_multiplicator = timeunit_day
    CASE ('hours')
      time_multiplicator = timeunit_hour
    CASE ('minuets')
      time_multiplicator = timeunit_minute
    CASE ('seconds')
      time_multiplicator = timeunit_second
    CASE DEFAULT
      STOP
    END SELECT

    DO i = 1, SIZE(time1)
      s0 = time1(i)*(1.0e-3_wp*time_multiplicator%ms)    ! in seconds
      s = s0
      days = 0.0_wp
      DO WHILE (s >= 86400.0_wp)
        s = s - 86400.0_wp
        days = days + 1.0_wp
      END DO
      DO WHILE (s < 0.0_wp)
        s = s + 86400.0_wp
        days = days - 1.0_wp
      END DO
      WRITE (0, '(i4,4f16.2)') i, time1(i), days, s, s/86400.0_wp
      s = s0
      days = FLOOR(s/86400.0_wp)
      ms = s - days*86400.0_wp
      WRITE (0, '(20x,3f16.2)') days, ms, ms/86400.0_wp
    END DO
  END BLOCK iteration1

  iteration2: BLOCK

    TYPE(juliandelta) :: time_multiplicator

    REAL(wp) :: s, s0, days, ms
    INTEGER :: i

    SELECT CASE (time2_unit)
    CASE ('days')
      time_multiplicator = timeunit_day
    CASE ('hours')
      time_multiplicator = timeunit_hour
    CASE ('minuets')
      time_multiplicator = timeunit_minute
    CASE ('seconds')
      time_multiplicator = timeunit_second
    CASE DEFAULT
      STOP
    END SELECT

    DO i = 1, SIZE(time2)
      s0 = time2(i)*(1.0e-3_wp*time_multiplicator%ms)    ! in seconds
      s = s0
      days = 0.0_wp
      DO WHILE (s >= 86400.0_wp)
        s = s - 86400.0_wp
        days = days + 1.0_wp
      END DO
      DO WHILE (s < 0.0_wp)
        s = s + 86400.0_wp
        days = days - 1.0_wp
      END DO
      WRITE (0, '(i4,4f16.2)') i, time2(i), days, s, s/86400.0_wp
      s = s0
      days = FLOOR(s/86400.0_wp)
      ms = s - days*86400.0_wp
      WRITE (0, '(20x,3f16.2)') days, ms, ms/86400.0_wp
    END DO
  END BLOCK iteration2

END PROGRAM test_cf_timeaxis
