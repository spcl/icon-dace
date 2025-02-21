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




MODULE mo_aes_convect_tables

  ! Lookup tables for convective adjustment code
  !
  ! D. Salmond, CRAY (UK), August 1991, original code
  ! L. Kornblueh, MPI, April 2003, cleanup and move of table setup code
  !                            from setphys.f90 in module
  ! M. Puetz, IBM, April 2009, added spline option
  ! L. Kornblueh, MPI, July 2010, corrected inconsistent use of saturation
  !                            vapour formulae
  ! L. Kornblueh, August 2010, corrected index bug causing inconsistencies
  !                            in table lookup for splines

  USE mo_kind,      ONLY: wp
  USE mo_exception, ONLY: message_text, message, finish
  USE mo_physical_constants, ONLY: alv, als, rd, rv, tmelt, cpd
  USE mo_aes_cop_config,     ONLY: aes_cop_config

  USE mo_model_domain,       ONLY: p_patch  ! for debugging only
  USE mo_math_constants,     ONLY: rad2deg  ! for debugging only

  IMPLICIT NONE

  SAVE

  PRIVATE

  ! Reference: Sonntag D., 1990: Important new values of the physical
  ! constants of 1986, vapour pressure formulations based on ITS-90,
  ! and psychrometer formulae. Z. Meteor. 70, pp 340-344.

  REAL(wp), PARAMETER :: cavl1 = -6096.9385_wp
  REAL(wp), PARAMETER :: cavl2 =    21.2409642_wp
  REAL(wp), PARAMETER :: cavl3 =    -2.711193_wp
  REAL(wp), PARAMETER :: cavl4 =     1.673952_wp
  REAL(wp), PARAMETER :: cavl5 =     2.433502_wp

  REAL(wp), PARAMETER :: cavi1 = -6024.5282_wp
  REAL(wp), PARAMETER :: cavi2 =    29.32707_wp
  REAL(wp), PARAMETER :: cavi3 =     1.0613868_wp
  REAL(wp), PARAMETER :: cavi4 =    -1.3198825_wp
  REAL(wp), PARAMETER :: cavi5 =    -0.49382577_wp

  ! variables public

  PUBLIC :: jptlucu1          ! lookup table lower bound
  PUBLIC :: jptlucu2          ! lookup table upper bound

  PUBLIC :: fjptlucu1         ! lookup table lower bound (real)
  PUBLIC :: fjptlucu2         ! lookup table upper bound (real)

  PUBLIC :: tlucua            ! table -- e_s*Rd/Rv
  PUBLIC :: tlucub            ! table -- derivative: d es/ d t
  PUBLIC :: tlucuc            ! table -- lv/cp

  PUBLIC :: tlucuaw           ! table -- e_s*Rd/Rv - forcing saturation wrt water
  Public :: tlucubw           ! table -- derivative: d es/ d t - forcing saturation wrt water
  PUBLIC :: tlucucw           ! table -- lv/cp - forcing saturation wrt water

  PUBLIC :: tlucu             ! fused table
  PUBLIC :: tlucuw            ! fused table

  PUBLIC :: lookupoverflow
  PUBLIC :: flookupoverflow

  PUBLIC :: prepare_ua_index_spline
  PUBLIC :: prepare_ua_index_spline_batch

  PUBLIC :: lookup_ua_spline
  PUBLIC :: lookup_ua_spline_batch

  PUBLIC :: lookup_uaw_spline
  PUBLIC :: lookup_ua_eor_uaw_spline
  PUBLIC :: lookup_ua_eor_uaw_spline_batch

  PUBLIC :: lookup_ua_list_spline
  PUBLIC :: lookup_uaw_list_spline

  PUBLIC :: prepare_ua_index

  PUBLIC :: lookup_ua
  PUBLIC :: lookup_uaw
  PUBLIC :: lookup_ua_eor_uaw

  PUBLIC :: lookup_ua_list
  PUBLIC :: lookup_uaw_list

  PUBLIC :: lookup_ubc
  PUBLIC :: lookup_ubc_list

  ! subroutines public

  PUBLIC :: init_convect_tables ! initialize LUTs
  PUBLIC :: lookuperror         ! error handling routine

  !---------------------------------------------------------------------------
  !
  ! Configuration section:
  !
  !-- temperature evaluation bounds:
  !
  REAL(wp), PARAMETER :: tlbound =  50.0_wp  ! lower bound [K]
  REAL(wp), PARAMETER :: tubound = 400.0_wp  ! upper bound [K]
  !
  !-- full table:
  !
  REAL(wp), PARAMETER :: fdeltat  =   0.001_wp  ! distance of table knots [K]
  ! division is not sufficiently precise, have to replace 1.0_wp/fdeltat
  REAL(wp), PARAMETER :: rfdeltat = 1000.0_wp
  !
  !-- spline interpolation:
  !
  REAL(wp), PARAMETER :: sdeltat  =   0.025_wp  ! distance of table knots [K]
  ! division is not sufficiently precise, have to replace 1.0_wp/sdeltat
  REAL(wp), PARAMETER :: rsdeltat =  40.0_wp
  !
  !---------------------------------------------------------------------------
  !
  ! Derived bounds and deltas:
  !
  !-- full table:
  !
  INTEGER,  PARAMETER :: jptlucu0 = NINT(rfdeltat*tmelt)   ! lookup table physical center
  INTEGER,  PARAMETER :: jptlucu1 = NINT(rfdeltat*tlbound) ! lookup table lower bound
  INTEGER,  PARAMETER :: jptlucu2 = NINT(rfdeltat*tubound) ! lookup table upper bound
  !
  !-- spline interpolation:
  !
  INTEGER,  PARAMETER :: lucupctr = NINT(rsdeltat*tmelt)   ! lookup table physical center
  INTEGER,  PARAMETER :: lucupmin = NINT(rsdeltat*tlbound) ! lookup table lower bound
  INTEGER,  PARAMETER :: lucupmax = NINT(rsdeltat*tubound) ! lookup table upper bound
  !
  !   For the spline interpolation the melting point tmelt = 273.15
  !   is resembled twice at index = lucubctr and index = lucubctr-1
  !   to be able to interpolate values close to the melting point
  !   correctly.
  !
  INTEGER, PARAMETER :: idxctr   = lucupctr
  INTEGER, PARAMETER :: idxctrm1 = idxctr-1
  !
  !---------------------------------------------------------------------------

  ! add float point constants for accelerated boundary guard code
  ! original big table
  REAL(wp), PARAMETER :: fjptlucu1 = REAL(jptlucu1,wp)
  REAL(wp), PARAMETER :: fjptlucu2 = REAL(jptlucu2,wp)
  ! new basic table bounds
  REAL(wp), PARAMETER :: flucupmin = REAL(lucupmin,wp)
  REAL(wp), PARAMETER :: flucupmax = REAL(lucupmax,wp)

  ! logical is fine and logical, but fp faster
  LOGICAL  :: lookupoverflow = .FALSE.          ! preset with false
  REAL(wp) :: flookupoverflow = 0.0_wp          ! preset with false

  REAL(wp) :: tlucua(jptlucu1-1:jptlucu2+1)     ! table - Es*Rd/Rv, mixed phases
  REAL(wp) :: tlucuaw(jptlucu1-1:jptlucu2+1)    ! table - Es*Rd/Rv, water phase only
  REAL(wp) :: tlucuad(jptlucu1-1:jptlucu2+1)    ! table - dEs/dT*Rd/Rv, mixed phases
  REAL(wp) :: tlucuawd(jptlucu1-1:jptlucu2+1)   ! table - dEs/dT*Rd/Rv, water phase only
  REAL(wp) :: tlucub(jptlucu1-1:jptlucu2+1)     ! table - inner dEs/dT*L/cp, mixed phases
  REAL(wp) :: tlucubw(jptlucu1-1:jptlucu2+1)    ! table - inner dEs/dT*L/cp, water phase only
  REAL(wp) :: tlucuc(jptlucu1-1:jptlucu2+1)     ! table - L/cp, mixed phases
  REAL(wp) :: tlucucw(jptlucu1-1:jptlucu2+1)    ! table - L/cp, water phase only
  !$ACC DECLARE CREATE(tlucua)

  ! fused tables for splines
  REAL(wp) :: tlucu(1:2,lucupmin-2:lucupmax+1)     ! fused table
  REAL(wp) :: tlucuw(1:2,lucupmin-2:lucupmax+1)    ! fused table

  !----------------------------------------------------------------------------
CONTAINS
  !----------------------------------------------------------------------------
  SUBROUTINE init_convect_tables
    REAL(wp) :: zalvdcp, zalsdcp
    REAL(wp) :: ztt, zldcp
    REAL(wp) :: zavm1, zavm2, zavm3, zavm4, zavm5
    REAL(wp) :: zminner,zdminner,zlinner,zdlinner
    INTEGER :: it

    ! check if tmelt is hit by one of the table entries (tricky to make right!):
    IF (ABS(REAL(lucupctr,wp)-rsdeltat*tmelt) > EPSILON(1.0_wp)) THEN
      CALL finish('init_convect_tables','distance of spline knots incorrect.')
    ENDIF

    WRITE (message_text,'(a,f8.2,a)')  'Melting temperature  : ', tmelt, ' [K]'
    CALL message('',message_text)
    WRITE (message_text,'(a,f8.2,a)')  'Lower table bound    : ', tlbound, ' [K]'
    CALL message('',message_text)
    WRITE (message_text,'(a,f8.2,a)')  'Upper table bound    : ', tubound, ' [K]'
    CALL message('',message_text)
    WRITE (message_text,'(a)')         'Full table settings  : '
    CALL message('',message_text)
    WRITE (message_text,'(a,f12.3,a)') '  Temperature delta  : ', fdeltat, ' [K]'
    CALL message('',message_text)
    WRITE (message_text,'(a,f12.3,a)') '  Inverse delta      : ', rfdeltat, ' [1/K]'
    CALL message('',message_text)
    WRITE (message_text,'(a,3i8)')     '  Indices [lb,ctr,ub]: ', jptlucu1, jptlucu0, jptlucu2
    CALL message('',message_text)
    WRITE (message_text,'(a)')         'Spline table settings: '
    CALL message('',message_text)
    WRITE (message_text,'(a,f12.3,a)') '  Temperature delta  : ', sdeltat, ' [K]'
    CALL message('',message_text)
    WRITE (message_text,'(a,f12.3,a)') '  Inverse delta      : ', rsdeltat, ' [1/K]'
    CALL message('',message_text)
    WRITE (message_text,'(a,3i8)')     '  Indices [lb,ctr,ub]: ', lucupmin, lucupctr, lucupmax
    CALL message('',message_text)
    WRITE (message_text,'(a,2i8)')     '  Table split indices: ', idxctr, idxctrm1
    CALL message('',message_text)

    zalvdcp = alv/cpd
    zalsdcp = als/cpd

    ! extend the table by one to safe guard code for accessing it+1, when
    ! it is already checked

    ! This are the old tables, extending over 350000 values: kills performance

!IBM* NOVECTOR
    DO it = jptlucu1-1, jptlucu2+1

      ztt = fdeltat*it
      IF ((ztt-tmelt) > 0.0_wp) THEN
        zldcp = zalvdcp
        zavm1 = cavl1
        zavm2 = cavl2
        zavm3 = cavl3
        zavm4 = cavl4
        zavm5 = cavl5
      ELSE
        zldcp = zalsdcp
        zavm1 = cavi1
        zavm2 = cavi2
        zavm3 = cavi3
        zavm4 = cavi4
        zavm5 = cavi5
      END IF

      zminner  = (zavm1/ztt+zavm2+zavm3*0.01_wp*ztt+zavm4*ztt*ztt*1.e-5_wp+zavm5*LOG(ztt))
      zdminner = (-zavm1/(ztt*ztt) + zavm3*0.01_wp + zavm4*ztt*2.e-5_wp + zavm5/ztt)

      zlinner  = (cavl1/ztt+cavl2+cavl3*0.01_wp*ztt+cavl4*ztt*ztt*1.e-5_wp+cavl5*LOG(ztt))
      zdlinner = (-cavl1/(ztt*ztt) + cavl3*0.01_wp + cavl4*ztt*2.e-5_wp + cavl5/ztt)

      tlucua(it)  = EXP(zminner)*rd/rv
      tlucuaw(it) = EXP(zlinner)*rd/rv

      tlucuad(it)  = zdminner*EXP(zminner)*rd/rv
      tlucuawd(it) = zdlinner*EXP(zlinner)*rd/rv

      tlucub(it)  = zldcp*zdminner
      tlucubw(it) = zalvdcp*zdlinner

      tlucuc(it)  = zldcp
      tlucucw(it) = zalvdcp

    END DO

    ! This are the new tables, base for the spline interpolation

!IBM* NOVECTOR
    DO it = lucupmin-2, idxctrm1
      ztt = sdeltat*(it+1)
      zldcp = zalsdcp
      zavm1 = cavi1
      zavm2 = cavi2
      zavm3 = cavi3
      zavm4 = cavi4
      zavm5 = cavi5

      zminner  = (zavm1/ztt+zavm2+zavm3*0.01_wp*ztt+zavm4*ztt*ztt*1.e-5_wp+zavm5*LOG(ztt))
      zdminner = (-zavm1/(ztt*ztt) + zavm3*0.01_wp + zavm4*ztt*2.e-5_wp + zavm5/ztt)

      zlinner  = (cavl1/ztt+cavl2+cavl3*0.01_wp*ztt+cavl4*ztt*ztt*1.e-5_wp+cavl5*LOG(ztt))
      zdlinner = (-cavl1/(ztt*ztt) + cavl3*0.01_wp + cavl4*ztt*2.e-5_wp + cavl5/ztt)

      tlucu(1,it)  = EXP(zminner)*rd/rv
      tlucu(2,it)  = sdeltat*zdminner*EXP(zminner)*rd/rv

      tlucuw(1,it) = EXP(zlinner)*rd/rv
      tlucuw(2,it) = sdeltat*zdlinner*EXP(zlinner)*rd/rv

    END DO

!IBM* NOVECTOR
    DO it = idxctr, lucupmax+1
      ztt = sdeltat*it
      zldcp = zalvdcp
      zavm1 = cavl1
      zavm2 = cavl2
      zavm3 = cavl3
      zavm4 = cavl4
      zavm5 = cavl5

      zminner  = (zavm1/ztt+zavm2+zavm3*0.01_wp*ztt+zavm4*ztt*ztt*1.e-5_wp+zavm5*LOG(ztt))
      zdminner = (-zavm1/(ztt*ztt) + zavm3*0.01_wp + zavm4*ztt*2.e-5_wp + zavm5/ztt)

      zlinner  = (cavl1/ztt+cavl2+cavl3*0.01_wp*ztt+cavl4*ztt*ztt*1.e-5_wp+cavl5*LOG(ztt))
      zdlinner = (-cavl1/(ztt*ztt) + cavl3*0.01_wp + cavl4*ztt*2.e-5_wp + cavl5/ztt)

      tlucu(1,it)  = EXP(zminner)*rd/rv
      tlucu(2,it)  = sdeltat*zdminner*EXP(zminner)*rd/rv

      tlucuw(1,it) = EXP(zlinner)*rd/rv
      tlucuw(2,it) = sdeltat*zdlinner*EXP(zlinner)*rd/rv
    END DO

    !$ACC ENTER DATA COPYIN(tlucu, tlucuw)
    !$ACC UPDATE DEVICE(tlucua) ASYNC(1)

  END SUBROUTINE init_convect_tables
  !----------------------------------------------------------------------------
  SUBROUTINE lookup_ubc(jcs,size, temp, ub, uc)
    INTEGER,            INTENT(in)  :: jcs, size
    REAL(wp),           INTENT(in)  :: temp(size)
    REAL(wp),           INTENT(out) :: ub(size)
    REAL(wp), OPTIONAL, INTENT(out) :: uc(size)

    REAL(wp) :: zalvdcp, zalsdcp, zavm1, zavm3, zavm4, zavm5, zldcp
    INTEGER :: jl

    zalvdcp = alv/cpd
    zalsdcp = als/cpd

    !$ACC DATA PRESENT(temp, ub)

    IF (PRESENT(uc)) THEN
      !$ACC DATA PRESENT(uc)
!IBM* NOVECTOR
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR PRIVATE(zavm1, zavm3, zavm4, zavm5, zldcp)
      DO jl = jcs, size
        zavm1 = MERGE(cavi1,cavl1,(tmelt-temp(jl)).GE.0._wp)
        zavm3 = MERGE(cavi3,cavl3,(tmelt-temp(jl)).GE.0._wp)
        zavm4 = MERGE(cavi4,cavl4,(tmelt-temp(jl)).GE.0._wp)
        zavm5 = MERGE(cavi5,cavl5,(tmelt-temp(jl)).GE.0._wp)
        zldcp = MERGE(zalsdcp,zalvdcp,(tmelt-temp(jl)).GE.0._wp)
        ub(jl) = zldcp*(-zavm1/(temp(jl)*temp(jl))+zavm3*0.01_wp+zavm4*temp(jl)*2.e-5_wp+zavm5/temp(jl))
        uc(jl) = zldcp
      END DO
      !$ACC END PARALLEL
      !$ACC END DATA
    ELSE
!IBM* NOVECTOR
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR PRIVATE(zavm1, zavm3, zavm4, zavm5, zldcp)
      DO jl = jcs, size
        zavm1 = MERGE(cavi1,cavl1,(tmelt-temp(jl)).GE.0._wp)
        zavm3 = MERGE(cavi3,cavl3,(tmelt-temp(jl)).GE.0._wp)
        zavm4 = MERGE(cavi4,cavl4,(tmelt-temp(jl)).GE.0._wp)
        zavm5 = MERGE(cavi5,cavl5,(tmelt-temp(jl)).GE.0._wp)
        zldcp = MERGE(zalsdcp,zalvdcp,(tmelt-temp(jl)).GE.0._wp)
        ub(jl) = zldcp*(-zavm1/(temp(jl)*temp(jl))+zavm3*0.01_wp+zavm4*temp(jl)*2.e-5_wp+zavm5/temp(jl))
      END DO
      !$ACC END PARALLEL
    END IF

    !$ACC END DATA

  END SUBROUTINE lookup_ubc
  !----------------------------------------------------------------------------
  SUBROUTINE lookup_ubc_list(jcs, size, kidx, list, temp, ub, uc)
    INTEGER,            INTENT(in)  :: jcs, size, kidx
    INTEGER,            INTENT(in)  :: list(kidx)
    REAL(wp),           INTENT(in)  :: temp(size)
    REAL(wp),           INTENT(out) :: ub(size)
    REAL(wp), OPTIONAL, INTENT(out) :: uc(size)

    REAL(wp) :: zalvdcp, zalsdcp, zavm1, zavm3, zavm4, zavm5, zldcp
    INTEGER :: nl, jl

    zalvdcp = alv/cpd
    zalsdcp = als/cpd

    !$ACC DATA PRESENT(list, temp, ub)

    IF (PRESENT(uc)) THEN
      !$ACC DATA PRESENT(uc)
!IBM* NOVECTOR
!IBM* ASSERT(NODEPS)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR PRIVATE(jl, zavm1, zavm3, zavm4, zavm5, zldcp)
      DO nl = jcs, kidx
        jl = list(nl)
        zavm1 = MERGE(cavi1,cavl1,(tmelt-temp(jl)).GE.0._wp)
        zavm3 = MERGE(cavi3,cavl3,(tmelt-temp(jl)).GE.0._wp)
        zavm4 = MERGE(cavi4,cavl4,(tmelt-temp(jl)).GE.0._wp)
        zavm5 = MERGE(cavi5,cavl5,(tmelt-temp(jl)).GE.0._wp)
        zldcp = MERGE(zalsdcp,zalvdcp,(tmelt-temp(jl)).GE.0._wp)
        ub(nl) = zldcp*(-zavm1/(temp(jl)*temp(jl))+zavm3*0.01_wp+zavm4*temp(jl)*2.e-5_wp+zavm5/temp(jl))
        uc(nl) = zldcp
      END DO
      !$ACC END PARALLEL
      !$ACC END DATA
    ELSE
!IBM* NOVECTOR
!IBM* ASSERT(NODEPS)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR PRIVATE(jl, zavm1, zavm3, zavm4, zavm5, zldcp)
      DO nl = jcs, kidx
        jl = list(nl)
        zavm1 = MERGE(cavi1,cavl1,(tmelt-temp(jl)).GE.0._wp)
        zavm3 = MERGE(cavi3,cavl3,(tmelt-temp(jl)).GE.0._wp)
        zavm4 = MERGE(cavi4,cavl4,(tmelt-temp(jl)).GE.0._wp)
        zavm5 = MERGE(cavi5,cavl5,(tmelt-temp(jl)).GE.0._wp)
        zldcp = MERGE(zalsdcp,zalvdcp,(tmelt-temp(jl)).GE.0._wp)
        ub(nl) = zldcp*(-zavm1/(temp(jl)*temp(jl))+zavm3*0.01_wp+zavm4*temp(jl)*2.e-5_wp+zavm5/temp(jl))
      END DO
      !$ACC END PARALLEL
    END IF

    !$ACC END DATA

  END SUBROUTINE lookup_ubc_list
  !----------------------------------------------------------------------------
  SUBROUTINE fetch_ua_spline(jcs,size,idx,zalpha,table,ua,dua)
    INTEGER,            INTENT(in)  :: jcs, size
    INTEGER,            INTENT(in)  :: idx(:)
    REAL(wp),           INTENT(in)  :: zalpha(:)
    REAL(wp),           INTENT(in)  :: table(1:2,lucupmin-2:lucupmax+1)
    REAL(wp), OPTIONAL, INTENT(out) :: ua(:), dua(:)

    REAL(wp) :: a, b, c, d, dx, ddx, x, bxa
    INTEGER :: jl

    !$ACC DATA PRESENT(idx, zalpha, table)

    IF (PRESENT(ua) .AND. .NOT. PRESENT(dua)) THEN
      !$ACC DATA PRESENT(ua)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR PRIVATE(x, dx, ddx, a, b, c, d, bxa)
      DO jl = jcs,size
        x = zalpha(jl)
        ! derivative and second derivative approximations (2 flops)
        dx   = table(1,idx(jl)+1) - table(1,idx(jl))
        ddx  = table(2,idx(jl)+1) + table(2,idx(jl))
        ! determine coefficients (2 fma + 1 flop)
        a = ddx - 2.0_wp*dx
        b = 3.0_wp*dx - ddx - table(2,idx(jl))
        c = table(2,idx(jl))
        d = table(1,idx(jl))
        ! Horner's scheme to compute the spline functions (3 fmas)
        bxa = b + x*a
        ua(jl) = d + x*(c + x*bxa)
      END DO
      !$ACC END PARALLEL
      !$ACC END DATA
    ELSE IF (PRESENT(ua) .AND. PRESENT(dua)) THEN
      !$ACC DATA PRESENT(ua, dua)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR PRIVATE(x, dx, ddx, a, b, c, d, bxa)
      DO jl = jcs,size
        x = zalpha(jl)
        ! derivate and second derivate approximations (2 flops)
        dx   = table(1,idx(jl)+1) - table(1,idx(jl))
        ddx  = table(2,idx(jl)+1) + table(2,idx(jl))
        ! determine coefficients (2 fma + 1 flop)
        a = ddx - 2.0_wp*dx
        b = 3.0_wp*dx - ddx - table(2,idx(jl))
        c = table(2,idx(jl))
        d = table(1,idx(jl))
        ! Horner's scheme to compute the spline functions (5 fmas + 1 flop)
        bxa = b + x*a
        ua(jl)  = d + x*(c + x*bxa)
        dua(jl) = rsdeltat*(c + x*(3.0_wp*bxa - b))
      END DO
      !$ACC END PARALLEL
      !$ACC END DATA
    END IF

    !$ACC END DATA

  END SUBROUTINE fetch_ua_spline
  !----------------------------------------------------------------------------
  SUBROUTINE fetch_ua(size, idx, table, u)
    INTEGER,  INTENT(in)  :: size
    INTEGER,  INTENT(in)  :: idx(size)
    REAL(wp), INTENT(in)  :: table(jptlucu1-1:jptlucu2+1)
    REAL(wp), INTENT(out) :: u(size)

    INTEGER :: jl

    DO jl = 1, size
      u(jl) = table(idx(jl))
    END DO

  END SUBROUTINE fetch_ua
  !----------------------------------------------------------------------------
  SUBROUTINE fetch_ua_list_spline(size, kidx1, kidx2, store_idx, lookup_idx, zalpha, table, ua, dua)
    INTEGER,            INTENT(in)    :: size, kidx1, kidx2
    INTEGER,            INTENT(in)    :: store_idx(size), lookup_idx(size)
    REAL(wp),           INTENT(in)    :: zalpha(size)
    REAL(wp),           INTENT(in)    :: table(1:2,lucupmin-2:lucupmax+1)
    REAL(wp), OPTIONAL, INTENT(inout) :: ua(size), dua(size)

    REAL(wp) :: a, b, c, d, dx, ddx, x, bxa
    INTEGER  :: jl, nl

    IF (PRESENT(ua) .AND. .NOT. PRESENT(dua)) THEN
      !$ACC DATA PRESENT(store_idx, lookup_idx, zalpha, table, ua)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR PRIVATE(jl, x, dx, ddx, a, b, c, d, bxa)
!IBM* ASSERT(NODEPS)
      DO nl = kidx1, kidx2
        jl = store_idx(nl)
        x = zalpha(jl)
        ! derivate and second derivate approximations (2 flops)
        dx   = table(1,lookup_idx(jl)+1) - table(1,lookup_idx(jl))
        ddx  = table(2,lookup_idx(jl)+1) + table(2,lookup_idx(jl))
        ! determine coefficients (2 fma + 1 flop)
        a = ddx - 2.0_wp*dx
        b = 3.0_wp*dx - ddx - table(2,lookup_idx(jl))
        c = table(2,lookup_idx(jl))
        d = table(1,lookup_idx(jl))
        ! Horner's scheme to compute the spline functions (3 fmas)
        bxa = b + x*a
        ua(jl) = d + x*(c + x*bxa)
      END DO
      !$ACC END PARALLEL
      !$ACC END DATA

    ELSE IF (PRESENT(ua) .AND. PRESENT(dua)) THEN
      !$ACC DATA PRESENT(store_idx, lookup_idx, zalpha, table, ua, dua)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR PRIVATE(jl, x, dx, ddx, a, b, c, d, bxa)
!IBM* ASSERT(NODEPS)
      DO nl = kidx1, kidx2
        jl = store_idx(nl)
        x = zalpha(jl)
        ! derivate and second derivate approximations (2 flops)
        dx   = table(1,lookup_idx(jl)+1) - table(1,lookup_idx(jl))
        ddx  = table(2,lookup_idx(jl)+1) + table(2,lookup_idx(jl))
        ! determine coefficients (2 fma + 1 flop)
        a = ddx - 2.0_wp*dx
        b = 3.0_wp*dx - ddx - table(2,lookup_idx(jl))
        c = table(2,lookup_idx(jl))
        d = table(1,lookup_idx(jl))
        ! Horner's scheme to compute the spline functions (5 fmas + 1 flop)
        bxa = b + x*a
        ua(jl)  = d + x*(c + x*bxa)
        dua(jl) = rsdeltat*(c + x*(3.0_wp*bxa - b))
      END DO
      !$ACC END PARALLEL
      !$ACC END DATA
    END IF

  END SUBROUTINE fetch_ua_list_spline
  !----------------------------------------------------------------------------
  SUBROUTINE ua_spline(table, x, idx, ua, dua)
    !$ACC ROUTINE SEQ
    REAL(wp),  INTENT(in)  :: table(1:2,lucupmin-2:lucupmax+1)
    REAL(wp),  INTENT(in)  :: x
    INTEGER,   INTENT(in)  :: idx
    REAL(wp),  INTENT(out) :: ua, dua

    REAL(wp) :: a, b, c, d, dx, ddx, bxa

    ! derivate and second derivate approximations (2 flops)
    dx   = table(1,idx+1) - table(1,idx)
    ddx  = table(2,idx+1) + table(2,idx)
    ! determine coefficients (2 fma + 1 flop)
    a = ddx - 2.0_wp*dx
    b = 3.0_wp*dx - ddx - table(2,idx)
    c = table(2,idx)
    d = table(1,idx)
    ! Horner's scheme to compute the spline functions (5 fmas + 1 flop)
    bxa = b + x*a
    ua  = d + x*(c + x*bxa)
    dua = rsdeltat*(c + x*(3.0_wp*bxa - b))
  END SUBROUTINE ua_spline

  SUBROUTINE fetch_ua_spline_batch(jcs,jce,batch_size,idx,zalpha,table,ua,dua)
    INTEGER,            INTENT(in)  :: jcs, jce, batch_size
    INTEGER,            INTENT(in)  :: idx(:, :)
    REAL(wp),           INTENT(in)  :: zalpha(:, :)
    REAL(wp),           INTENT(in)  :: table(1:2,lucupmin-2:lucupmax+1)
    REAL(wp), OPTIONAL, INTENT(out) :: ua(:, :), dua(:, :)

    REAL(wp) :: mydua
    INTEGER :: jl, batch
    LOGICAL :: need_dua

    need_dua = PRESENT(dua)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) NO_CREATE(dua) GANG VECTOR COLLAPSE(2) ASYNC(1)
    DO batch = 1,batch_size
      DO jl = jcs,jce
        CALL ua_spline(table, zalpha(jl,batch), idx(jl,batch), ua(jl,batch), mydua)
        IF ( need_dua ) dua(jl,batch) = mydua
      END DO
    END DO

  END SUBROUTINE fetch_ua_spline_batch

  SUBROUTINE fetch_ua_list_spline_batch(jcs, jce, batch_size,          &
                                        lookup_idx, iphase,            &
                                        zalpha, table, tablew, ua, dua)

    INTEGER,            INTENT(in)    :: jcs,jce,batch_size
    INTEGER,            INTENT(in)    :: lookup_idx(:,:), iphase(:,:)
    REAL(wp),           INTENT(in)    :: zalpha(:,:)
    REAL(wp),           INTENT(in)    :: table(1:2,lucupmin-2:lucupmax+1)

    ! For merging cases with temp >= t_melt and temp < t_melt
    REAL(wp),           INTENT(in)    :: tablew(1:2,lucupmin-2:lucupmax+1)
    REAL(wp), OPTIONAL, INTENT(inout) :: ua(:,:), dua(:,:)

    INTEGER  :: jl, batch, idx
    LOGICAL  :: need_dua
    REAL(wp) :: myua, mydua, x

    need_dua = PRESENT(dua)

    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) NO_CREATE(dua) COLLAPSE(2) ASYNC(1)
    DO batch = 1, batch_size
      DO jl = jcs, jce

        idx = lookup_idx(jl,batch)
        x   = zalpha    (jl,batch)
        IF (iphase(jl,batch) == 1) THEN
          CALL ua_spline(table,  x, idx, myua, mydua)
        ELSE
          CALL ua_spline(tablew, x, idx, myua, mydua)
        END IF

        ua(jl,batch) = myua
        IF ( need_dua ) dua(jl,batch) = mydua
      END DO
    END DO

  END SUBROUTINE fetch_ua_list_spline_batch
  !----------------------------------------------------------------------------
  SUBROUTINE fetch_ua_list(size, kidx, store_idx, lookup_idx, table, u)
    INTEGER,  INTENT(in)    :: size, kidx
    INTEGER,  INTENT(in)    :: store_idx(kidx), lookup_idx(size)
    REAL(wp), INTENT(in)    :: table(jptlucu1-1:jptlucu2+1)
    REAL(wp), INTENT(inout) :: u(size)

    INTEGER :: jl, nl

!IBM* ASSERT(NODEPS)
    DO nl = 1, kidx
      jl = store_idx(nl)
      u(jl) = table(lookup_idx(jl))
    END DO

  END SUBROUTINE fetch_ua_list
  !----------------------------------------------------------------------------
  SUBROUTINE lookup_ua_eor_uaw_spline(jcs, size, idx, zalpha, nphase, iphase, ua, dua)
    INTEGER,            INTENT(in)  :: jcs, size, nphase
    INTEGER,            INTENT(in)  :: idx(size), iphase(size)
    REAL(wp),           INTENT(in)  :: zalpha(size)
    REAL(wp), OPTIONAL, INTENT(out) :: ua(size), dua(size)

    INTEGER :: tmpidx(size)
    INTEGER :: jl,iw,inw

!    REAL(wp) :: ua1(size), dua1(size)

    IF (nphase == 0) THEN ! temp > 0 for all indices
      IF (PRESENT(dua)) THEN
        CALL fetch_ua_spline(jcs, size, idx, zalpha, tlucuw, ua, dua)
      ELSE
        CALL fetch_ua_spline(jcs, size, idx, zalpha, tlucuw, ua)
      END IF

    ELSE IF (nphase == size-jcs+1) THEN ! temp < 0 for all indices
      IF (PRESENT(dua)) THEN
        CALL fetch_ua_spline(jcs, size, idx, zalpha, tlucu, ua, dua)
      ELSE
        CALL fetch_ua_spline(jcs, size, idx, zalpha, tlucu, ua)
      END IF

    ELSE
      !$ACC DATA PRESENT(idx, iphase) &
      !$ACC   CREATE(tmpidx)
      ! mixed case, must build store index
      iw = jcs    ! store indices with temp < tmelt at iw (iw = jcs:jcs+nphase-1)
      inw = size  ! store indices with temp >= tmelt at inw (inw = jcs+nphase : size)
      !$ACC UPDATE HOST(iphase) ASYNC(1)
      !$ACC WAIT(1)
      DO jl = jcs,size
        tmpidx(iw)  = jl  ! lower part of tmpidx() filled with cond = .true.
        tmpidx(inw) = jl  ! upper part of tmpidx() filled with cond = .false.
        iw = iw + iphase(jl)          ! iphase(jl)=1 if temp(jl) < tmelt, =0 if temp .ge. tmelt
        inw = inw - (1 - iphase(jl))
      END DO
      !$ACC UPDATE DEVICE(tmpidx) ASYNC(1)
      iw = iw - 1

      IF (PRESENT(dua)) THEN
        CALL fetch_ua_list_spline(size, jcs, iw+1, tmpidx, idx, zalpha, tlucu , ua=ua, dua=dua)
        CALL fetch_ua_list_spline(size, iw+1, size, tmpidx, idx, zalpha, tlucuw, ua=ua, dua=dua)
      ELSE
        CALL fetch_ua_list_spline(size, jcs, iw+1, tmpidx, idx, zalpha, tlucu , ua=ua)
        CALL fetch_ua_list_spline(size, iw+1, size, tmpidx, idx, zalpha, tlucuw, ua=ua)
      END IF
      !$ACC END DATA
    END IF

  END SUBROUTINE lookup_ua_eor_uaw_spline
  !----------------------------------------------------------------------------
  SUBROUTINE lookup_ua_eor_uaw_spline_batch(jcs, jce, batch_size, idx, iphase, zalpha, ua, dua)
    INTEGER,            INTENT(in)    :: jcs, jce, batch_size
    INTEGER,            INTENT(in)    :: idx(:,:), iphase(:,:)
    REAL(wp),           INTENT(in)    :: zalpha(:,:)
    REAL(wp), OPTIONAL, INTENT(inout) :: ua(:,:), dua(:,:)

    IF (PRESENT(dua)) THEN
      CALL fetch_ua_list_spline_batch(jcs, jce, batch_size, idx, iphase, zalpha, tlucu, tlucuw, ua=ua, dua=dua)
    ELSE
      CALL fetch_ua_list_spline_batch(jcs, jce, batch_size, idx, iphase, zalpha, tlucu, tlucuw, ua=ua)
    END IF

  END SUBROUTINE lookup_ua_eor_uaw_spline_batch
  !----------------------------------------------------------------------------
  SUBROUTINE lookup_ua_eor_uaw(size, idx, nphase, iphase, ua, dua)
    INTEGER,            INTENT(in)    :: size, nphase
    INTEGER,            INTENT(in)    :: idx(size), iphase(size)
    REAL(wp),           INTENT(inout) :: ua(size)
    REAL(wp), OPTIONAL, INTENT(inout) :: dua(size)

    INTEGER :: tmpidx(size)
    INTEGER :: jl, iw, inw

    IF (nphase == 0) THEN
      ! only water
      CALL fetch_ua(size,idx,tlucuaw,ua)
      IF (PRESENT(dua)) CALL fetch_ua(size,idx,tlucuawd,dua)
    ELSE IF (nphase == size) THEN
      ! only ice
      CALL fetch_ua(size,idx,tlucua,ua)
      IF (PRESENT(dua)) CALL fetch_ua(size,idx,tlucuad,dua)
    ELSE
      ! mixed case, must build store index
      iw = 1
      inw = size
      DO jl = 1,size
        tmpidx(iw)  = jl  ! lower part of tmpidx() filled with cond = .true.
        tmpidx(inw) = jl  ! upper part of tmpidx() filled with cond = .false.
        iw = iw + iphase(jl)
        inw = inw - (1 - iphase(jl))
      END DO
      iw = iw - 1
      CALL fetch_ua_list(size, iw     , tmpidx(1)   , idx, tlucua, ua)
      CALL fetch_ua_list(size, size-iw, tmpidx(iw+1), idx, tlucuaw, ua)
      IF (PRESENT(dua)) THEN
        CALL fetch_ua_list(size, iw     , tmpidx(1)   , idx, tlucuad, dua)
        CALL fetch_ua_list(size, size-iw, tmpidx(iw+1), idx, tlucuawd, dua)
      END IF
    END IF

  END SUBROUTINE lookup_ua_eor_uaw
  !----------------------------------------------------------------------------
  SUBROUTINE lookup_ua_spline(jcs, size, idx, zalpha, ua, dua)
    INTEGER,            INTENT(in)  :: jcs, size
    INTEGER,            INTENT(in)  :: idx(size)
    REAL(wp),           INTENT(in)  :: zalpha(size)
    REAL(wp), OPTIONAL, INTENT(out) :: ua(size), dua(size)

    CALL fetch_ua_spline(jcs, size, idx, zalpha, tlucu, ua, dua)

  END SUBROUTINE lookup_ua_spline
  !----------------------------------------------------------------------------
  SUBROUTINE lookup_ua_spline_batch(jcs, jce, batch_size, idx, zalpha, ua, dua)
    INTEGER,            INTENT(in)  :: jcs, jce, batch_size
    INTEGER,            INTENT(in)  :: idx(:,:)
    REAL(wp),           INTENT(in)  :: zalpha(:,:)
    REAL(wp), OPTIONAL, INTENT(out) :: ua(:,:), dua(:,:)

    CALL fetch_ua_spline_batch(jcs, jce, batch_size, idx, zalpha, tlucu, ua, dua)

  END SUBROUTINE lookup_ua_spline_batch
  !----------------------------------------------------------------------------
  SUBROUTINE lookup_ua(size, idx, ua, dua)
    INTEGER,            INTENT(in)    :: size
    INTEGER,            INTENT(in)    :: idx(size)
    REAL(wp),           INTENT(inout) :: ua(size)
    REAL(wp), OPTIONAL, INTENT(inout) :: dua(size)

    CALL fetch_ua(size, idx, tlucua, ua)
    IF (PRESENT(dua)) CALL fetch_ua(size, idx, tlucuad, dua)

  END SUBROUTINE lookup_ua
  !----------------------------------------------------------------------------
  SUBROUTINE lookup_uaw_spline(jcs, size, idx, zalpha, ua, dua)
    INTEGER,            INTENT(in)  :: jcs, size
    INTEGER,            INTENT(in)  :: idx(size)
    REAL(wp),           INTENT(in)  :: zalpha(size)
    REAL(wp), OPTIONAL, INTENT(out) :: ua(size), dua(size)

    CALL fetch_ua_spline(jcs, size,idx,zalpha,tlucuw,ua,dua)

  END SUBROUTINE lookup_uaw_spline
  !----------------------------------------------------------------------------
  SUBROUTINE lookup_uaw(size, idx, ua, dua)
    INTEGER,            INTENT(in)  :: size
    INTEGER,            INTENT(in)  :: idx(size)
    REAL(wp), OPTIONAL, INTENT(out) :: ua(size), dua(size)

    CALL fetch_ua(size, idx, tlucuaw, ua)
    IF (PRESENT(dua)) CALL fetch_ua(size, idx, tlucuawd, dua)

  END SUBROUTINE lookup_uaw
  !----------------------------------------------------------------------------
SUBROUTINE prepare_ua_index_spline(jg, name, jcs, size, temp, idx, zalpha, &
    &                                xi, nphase, zphase, iphase,             &
    &                                klev, kblock, kblock_size)
    INTEGER,            INTENT(in)  :: jg
    CHARACTER(len=*),   INTENT(in)  :: name
    INTEGER,            INTENT(in)  :: size, jcs ! start and end index of block
    REAL(wp),           INTENT(in)  :: temp(size)
    INTEGER,            INTENT(out) :: idx(size)
    REAL(wp),           INTENT(out) :: zalpha(size)

    INTEGER,  OPTIONAL, INTENT(in)  :: klev
    INTEGER,  OPTIONAL, INTENT(in)  :: kblock
    INTEGER,  OPTIONAL, INTENT(in)  :: kblock_size

    REAL(wp), OPTIONAL, INTENT(in)  :: xi(size)
    INTEGER,  OPTIONAL, INTENT(out) :: nphase
    REAL(wp), OPTIONAL, INTENT(out) :: zphase(size)
    INTEGER,  OPTIONAL, INTENT(out) :: iphase(size)

    REAL(wp) :: ztt, ztshft, zinbounds, ztmin,ztmax,znphase,ztest
    INTEGER :: jl

    ! Shortcuts to components of aes_cop_config
    !
    REAL(wp) :: csecfrl, cthomi
    !
    !$ACC DATA PRESENT(temp, idx, zalpha) &
    !$ACC   NO_CREATE(xi, zphase, iphase)
    !
    csecfrl = aes_cop_config(jg)%csecfrl
    cthomi  = aes_cop_config(jg)%cthomi

    zinbounds = 1._wp
    ztmin = flucupmin
    ztmax = flucupmax
    IF (PRESENT(xi)) THEN
      znphase = 0.0_wp
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR PRIVATE(ztshft, ztt, ztest) REDUCTION(+: znphase) PRIVATE(zinbounds)
      DO jl = jcs,size
        ztshft = MERGE(1.0_wp,0.0_wp,(tmelt-temp(jl)).GE.0._wp)
        ztt = rsdeltat*temp(jl)
        zalpha(jl) = ztt - AINT(ztt)
        idx(jl) = INT(ztt-ztshft)
        zinbounds = MERGE(0.0_wp,zinbounds,(ztmin-ztt).GE.0._wp)
        zinbounds = MERGE(0.0_wp,zinbounds,(ztt-ztmax).GE.0._wp)

        ! check dual phase conditions:
        ! lo2 = (ptm1(jl,jk) < cthomi) .OR. (ptm1(jl,jk) < tmelt .AND. zxised > csecfrl)
        ztest = MERGE(0.0_wp,1.0_wp,(temp(jl)-tmelt).GE.0._wp)
        ztest = MERGE(0.0_wp,ztest,(csecfrl-xi(jl)).GE.0._wp)
        ztest = MERGE(ztest,1.0_wp,(temp(jl)-cthomi).GE.0._wp)
        ! normalize ztest to 0 and 1
        iphase(jl) = INT(ztest)
        zphase(jl) = ztest-0.5_wp
        znphase = znphase + ztest
      END DO
      !$ACC END PARALLEL
      nphase = INT(znphase)
    ELSE
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR PRIVATE(ztshft, ztt) PRIVATE(zinbounds)
      DO jl = jcs, size
        ztshft = MERGE(1.0_wp,0.0_wp,(tmelt-temp(jl)).GE.0._wp)
        ztt = rsdeltat*temp(jl)
        zalpha(jl) = ztt - AINT(ztt)
        idx(jl) = INT(ztt-ztshft)
        zinbounds = MERGE(0.0_wp,zinbounds,(ztmin-ztt).GE.0._wp)
        zinbounds = MERGE(0.0_wp,zinbounds,(ztt-ztmax).GE.0._wp)
      END DO
      !$ACC END PARALLEL
    END IF

    ! if one index was out of bounds -> print error and exit
    IF (zinbounds == 0.0_wp) THEN
      IF ( PRESENT(kblock) .AND. PRESENT(kblock_size) .AND. PRESENT(klev) ) THEN
        !$ACC UPDATE HOST(temp) ASYNC(1)
        !$ACC WAIT(1)
        ! tied to patch(1), does not yet work for nested grids
        DO jl = 1, size
          ztt = rsdeltat*temp(jl)
          IF ( ztt <= ztmin .OR. ztt >= ztmax ) THEN
            WRITE ( 0 , '(a,a,a,a,i5,a,i8,a,f8.2,a,f8.2,a,f8.2)' )                                 &
                 & ' Lookup table problem in ', TRIM(name), ' at ',                                &
                 & ' level   =',klev,                                                              &
                 & ' cell ID =',p_patch(1)%cells%decomp_info%glb_index((kblock-1)*kblock_size+jl), &
                 & ' lon(deg)=',p_patch(1)%cells%center(jl,kblock)%lon*rad2deg,                    &
                 & ' lat(deg)=',p_patch(1)%cells%center(jl,kblock)%lat*rad2deg,                    &
                 & ' value   =',temp(jl)
          ENDIF
        ENDDO
      ENDIF
      CALL lookuperror(name, 'prepare_ua_index_spline')
    END IF

    !$ACC WAIT(1)
    !$ACC END DATA

  END SUBROUTINE prepare_ua_index_spline
  !----------------------------------------------------------------------------
  SUBROUTINE prepare_ua_index_spline_batch(name, jcs, jce, batch_size,                &
    &                                      temp, idx, zalpha,                         &
    &                                      xi, nphase, zphase, iphase,                &
    &                                      kblock, kblock_size, opt_need_host_nphase, &
    &                                      opt_sanitize_index, csecfrl, cthomi)
    CHARACTER(len=*),   INTENT(in)    :: name
    INTEGER,            INTENT(in)    :: jcs, jce, batch_size ! start and end index of block
    REAL(wp),           INTENT(in)    :: temp(:,:)
    INTEGER,            INTENT(inout) :: idx(:,:)
    REAL(wp),           INTENT(inout) :: zalpha(:,:)

    INTEGER,  OPTIONAL, INTENT(in)    :: kblock
    INTEGER,  OPTIONAL, INTENT(in)    :: kblock_size

    REAL(wp), OPTIONAL, INTENT(in)    :: xi(:,:)
    INTEGER,  OPTIONAL, INTENT(inout) :: nphase(:)
    REAL(wp), OPTIONAL, INTENT(inout) :: zphase(:,:)
    INTEGER,  OPTIONAL, INTENT(inout) :: iphase(:,:)

    LOGICAL,  OPTIONAL, INTENT(in)  :: opt_need_host_nphase, opt_sanitize_index
    REAL(wp), INTENT(IN) :: csecfrl ! [kg/kg]  minimum in-cloud water mass mixing ratio in mixed phase clouds
    REAL(wp), INTENT(IN) :: cthomi ! [K]      maximum temperature for homogeneous freezing

    LOGICAL :: sanitize_index

    REAL(wp) :: ztt, ztshft, ztmin,ztmax,znphase(batch_size),ztest,myznphase
    INTEGER :: jl, jgang, jvec, batch, zoutofbounds, myzoutofbounds
    INTEGER :: zoutofbounds_vec(batch_size)
    INTEGER, PARAMETER :: tile_size = 1024

    !
    !$ACC DATA PRESENT(temp, idx, zalpha)
    !

    zoutofbounds = 0
    ztmin = flucupmin
    ztmax = flucupmax

    sanitize_index = .FALSE.
    IF (PRESENT(opt_sanitize_index)) THEN
      IF (opt_sanitize_index) THEN
        sanitize_index = .TRUE.
      END IF
    END IF

    IF (PRESENT(xi)) THEN
      !$ACC DATA PRESENT(xi, zphase, iphase, nphase) &
      !$ACC   CREATE(znphase, zoutofbounds_vec)

      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1)
      DO batch = 1,batch_size
        znphase(batch) = 0.0_wp
        zoutofbounds_vec(batch) = 0.0_wp
      END DO

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG COLLAPSE(2)
      DO batch = 1,batch_size
        DO jgang = jcs, jce, tile_size

          myznphase = 0.0_wp
          myzoutofbounds = 0
          !$ACC LOOP VECTOR REDUCTION(+: myznphase, myzoutofbounds)
          DO jvec = 0,tile_size-1
            jl = jgang + jvec
            IF (jl > jce) CYCLE

            ztshft = MERGE(1.0_wp,0.0_wp,(tmelt-temp(jl,batch)).GE.0._wp)
            ztt = rsdeltat*temp(jl,batch)
            zalpha(jl,batch) = ztt - AINT(ztt)
            idx(jl,batch) = INT(ztt-ztshft)
            IF ( ztt <= ztmin .OR. ztt >= ztmax ) myzoutofbounds = myzoutofbounds + 1

            ! check dual phase conditions:
            ! lo2 = (ptm1(jl,jk) < cthomi) .OR. (ptm1(jl,jk) < tmelt .AND. zxised > csecfrl)
            ztest = MERGE(0.0_wp,1.0_wp,(temp(jl,batch)-tmelt).GE.0._wp)
            ztest = MERGE(0.0_wp,ztest,(csecfrl-xi(jl,batch)).GE.0._wp)
            ztest = MERGE(ztest,1.0_wp,(temp(jl,batch)-cthomi).GE.0._wp)
            ! normalize ztest to 0 and 1
            iphase(jl,batch) = INT(ztest)
            zphase(jl,batch) = ztest-0.5_wp

            myznphase = myznphase + ztest
          END DO

          !$ACC ATOMIC
          znphase(batch) = znphase(batch) + myznphase

          IF (zoutofbounds_vec(batch) > 0) THEN
            !$ACC ATOMIC
            zoutofbounds_vec(batch) = zoutofbounds_vec(batch) + myzoutofbounds
          END IF

        END DO
      END DO
      !$ACC END PARALLEL

      IF (sanitize_index) THEN
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) REDUCTION(+: zoutofbounds) ASYNC(1)
        DO batch = 1,batch_size
          nphase = INT(znphase)
          zoutofbounds = zoutofbounds + zoutofbounds_vec(batch)
        END DO
      END IF

      IF ( PRESENT(opt_need_host_nphase) ) THEN
        IF ( opt_need_host_nphase ) THEN
          !$ACC UPDATE HOST(nphase) ASYNC(1)
          !$ACC WAIT(1)
        END IF
      END IF

      !$ACC END DATA

    ELSE

      IF (sanitize_index) THEN
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) REDUCTION(+: zoutofbounds) GANG VECTOR COLLAPSE(2) ASYNC(1)
        DO batch = 1,batch_size
          DO jl = jcs, jce
            ztshft = MERGE(1.0_wp,0.0_wp,(tmelt-temp(jl,batch)).GE.0._wp)
            ztt = rsdeltat*temp(jl,batch)
            zalpha(jl,batch) = ztt - AINT(ztt)
            idx(jl,batch) = INT(ztt-ztshft)
            IF ( ztt <= ztmin .OR. ztt >= ztmax ) zoutofbounds = zoutofbounds + 1
          END DO
        END DO
      ELSE
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
        DO batch = 1,batch_size
          DO jl = jcs, jce
            ztshft = MERGE(1.0_wp,0.0_wp,(tmelt-temp(jl,batch)).GE.0._wp)
            ztt = rsdeltat*temp(jl,batch)
            zalpha(jl,batch) = ztt - AINT(ztt)
            idx(jl,batch) = INT(ztt-ztshft)
          END DO
        END DO
      END IF ! sanitize

    END IF ! present(xi)

    ! if one index was out of bounds -> print error and exit
    IF (zoutofbounds > 0) THEN
      IF ( PRESENT(kblock) .AND. PRESENT(kblock_size) ) THEN
        !$ACC UPDATE HOST(temp) ASYNC(1)
        !$ACC WAIT(1)
        ! tied to patch(1), does not yet work for nested grids
        DO batch = 1,batch_size
          DO jl = jcs, jce
            ztt = rsdeltat*temp(jl,batch)
            IF ( ztt <= ztmin .OR. ztt >= ztmax ) THEN
              WRITE ( 0 , '(a,a,a,a,i5,a,i8,a,f8.2,a,f8.2,a,f8.2)' )                               &
                 & ' Lookup table problem in ', TRIM(name), ' at ',                                &
                 & ' level   =',batch,                                                             &
                 & ' cell ID =',p_patch(1)%cells%decomp_info%glb_index((kblock-1)*kblock_size+jl), &
                 & ' lon(deg)=',p_patch(1)%edges%center(jl,kblock)%lon*rad2deg,                    &
                 & ' lat(deg)=',p_patch(1)%edges%center(jl,kblock)%lat*rad2deg,                    &
                 & ' value   =',temp(jl,batch)
            ENDIF
          ENDDO
        ENDDO
      ENDIF
      CALL lookuperror(name, 'prepare_ua_index_spline_batch')
    END IF

    !$ACC END DATA

  END SUBROUTINE prepare_ua_index_spline_batch
  !----------------------------------------------------------------------------
  SUBROUTINE prepare_ua_index(jg, name, size, temp, idx, xi, nphase, zphase, iphase)
    INTEGER,            INTENT(in)  :: jg
    CHARACTER(len=*),   INTENT(in)  :: name
    INTEGER,            INTENT(in)  :: size
    REAL(wp),           INTENT(in)  :: temp(size)
    INTEGER,            INTENT(out) :: idx(size)
    REAL(wp), OPTIONAL, INTENT(in)  :: xi(size)
    INTEGER,  OPTIONAL, INTENT(out) :: nphase
    REAL(wp), OPTIONAL, INTENT(out) :: zphase(size)
    INTEGER,  OPTIONAL, INTENT(out) :: iphase(size)

    REAL(wp) :: ztt, zinbounds, ztmin,ztmax,znphase,ztest
    INTEGER :: jl

    ! Shortcuts to components of aes_cop_config
    !
    REAL(wp), POINTER :: csecfrl, cthomi
    !
    csecfrl => aes_cop_config(jg)% csecfrl
    cthomi  => aes_cop_config(jg)% cthomi

    ! first compute all lookup indices and check if they are all within allowed bounds

    zinbounds = 1.0_wp
    ztmin = fjptlucu1
    ztmax = fjptlucu2
    IF (PRESENT(xi)) THEN
      znphase = 0.0_wp
      DO jl = 1, size
        ztt = ANINT(rfdeltat*temp(jl))
        idx(jl) = INT(ztt)
        zinbounds = MERGE(0.0_wp,zinbounds,(ztmin-ztt).GE.0._wp)
        zinbounds = MERGE(0.0_wp,zinbounds,(ztt-ztmax).GE.0._wp)
        ! check dual phase conditions
        ! lo2 = (ptm1(jl,jk) < cthomi) .OR. (ptm1(jl,jk) < tmelt .AND. zxised > csecfrl)
        ztest = MERGE(0.0_wp,1.0_wp,(temp(jl)-tmelt).GE.0._wp)
        ztest = MERGE(0.0_wp,ztest,(csecfrl-xi(jl)).GE.0._wp)
        ztest = MERGE(ztest,1.0_wp,(temp(jl)-cthomi).GE.0._wp)
        ! normalize ztest to 0 and 1
        iphase(jl) = INT(ztest)
        zphase(jl) = ztest-0.5_wp
        znphase = znphase + ztest
      END DO
      nphase = INT(znphase)
    ELSE
      DO jl = 1,size
        ztt = ANINT(rfdeltat*temp(jl))
        idx(jl) = INT(ztt)
        zinbounds = MERGE(0.0_wp,zinbounds,(ztmin-ztt).GE.0._wp)
        zinbounds = MERGE(0.0_wp,zinbounds,(ztt-ztmax).GE.0._wp)
      END DO
    END IF
    ! if one index was out of bounds -> print error and exit
    IF (zinbounds == 0.0_wp) CALL lookuperror(name, 'prepare_ua_index')

  END SUBROUTINE prepare_ua_index
  !----------------------------------------------------------------------------
  SUBROUTINE lookup_ua_list_spline(name, jcs, size, kidx, list, temp, ua, dua, &
    &                              klev, kblock, kblock_size)
    CHARACTER(len=*),   INTENT(in)  :: name
    INTEGER,            INTENT(in)  :: jcs, size, kidx
    INTEGER,            INTENT(in)  :: list(:)
    REAL(wp),           INTENT(in)  :: temp(:)
    REAL(wp), OPTIONAL, INTENT(out) :: ua(:), dua(:)
    INTEGER,  OPTIONAL, INTENT(in)  :: klev
    INTEGER,  OPTIONAL, INTENT(in)  :: kblock
    INTEGER,  OPTIONAL, INTENT(in)  :: kblock_size

    INTEGER :: idx(size)
    REAL(wp) :: zalpha(size)

    REAL(wp) :: ztt, ztshft, zinbounds, ztmax, ztmin
    INTEGER :: nl, jl

    !$ACC DATA PRESENT(list, temp) &
    !$ACC   NO_CREATE(ua, dua) &
    !$ACC   CREATE(idx, zalpha)

    zinbounds = 1.0_wp
    ztmin = flucupmin
    ztmax = flucupmax

    ! first compute all lookup indices and check if they are all within allowed bounds

!IBM* ASSERT(NODEPS)
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR PRIVATE(jl, ztshft, ztt) PRIVATE(zinbounds)
    DO nl = jcs, kidx
      jl = list(nl)
      ztshft = MERGE(1.0_wp,0.0_wp,(tmelt-temp(jl)).GE.0._wp)
      ztt = rsdeltat*temp(jl)
      zalpha(nl) = ztt - AINT(ztt)
      idx(nl) = INT(ztt-ztshft)
      zinbounds = MERGE(0.0_wp,zinbounds,(ztmin-ztt).GE.0._wp)
      zinbounds = MERGE(0.0_wp,zinbounds,(ztt-ztmax).GE.0._wp)
    END DO
    !$ACC END PARALLEL

    ! if one index was out of bounds -> print error and exit
    IF (zinbounds == 0.0_wp) THEN
      IF ( PRESENT(kblock) .AND. PRESENT(kblock_size) .AND. PRESENT(klev) ) THEN
        ! tied to patch(1), does not yet work for nested grids
        !$ACC UPDATE HOST(temp) ASYNC(1)
        !$ACC WAIT(1)
        DO jl = 1, size
          ztt = rsdeltat*temp(jl)
          IF ( ztt <= ztmin .OR. ztt >= ztmax ) THEN
            WRITE ( 0 , '(a,a,a,a,i5,a,i8,a,f8.2,a,f8.2,a,f8.2)' )                                 &
                 & ' Lookup table problem in ', TRIM(name), ' at ',                                &
                 & ' level   =',klev,                                                              &
                 & ' cell ID =',p_patch(1)%cells%decomp_info%glb_index((kblock-1)*kblock_size+jl), &
                 & ' lon(deg)=',p_patch(1)%cells%center(jl,kblock)%lon*rad2deg,                    &
                 & ' lat(deg)=',p_patch(1)%cells%center(jl,kblock)%lat*rad2deg,                    &
                 & ' value   =',temp(jl)
          ENDIF
        ENDDO
      ENDIF
      CALL lookuperror(name, 'lookup_ua_list_spline')
    ENDIF
    CALL fetch_ua_spline(jcs, kidx, idx, zalpha, tlucu, ua, dua)

    !$ACC WAIT(1)
    !$ACC END DATA
  END SUBROUTINE lookup_ua_list_spline
  !----------------------------------------------------------------------------
  SUBROUTINE lookup_ua_list(name, size, kidx, list, temp, ua, dua)
    CHARACTER(len=*),   INTENT(in)  :: name
    INTEGER,            INTENT(in)  :: size, kidx
    INTEGER,            INTENT(in)  :: list(kidx)
    REAL(wp),           INTENT(in)  :: temp(size)
    REAL(wp), OPTIONAL, INTENT(out) :: ua(size), dua(size)

    INTEGER :: idx(size)

    REAL(wp) :: ztt, zinbounds, ztmax, ztmin
    INTEGER :: nl, jl

    zinbounds = 1.0_wp
    ztmin = fjptlucu1
    ztmax = fjptlucu2

    ! first compute all lookup indices and check if they are all within allowed bounds

!IBM* ASSERT(NODEPS)
    DO nl = 1,kidx
      jl = list(nl)
      ztt = ANINT(rfdeltat*temp(jl))
      idx(nl) = INT(ztt)
      zinbounds = MERGE(0.0_wp,zinbounds,(ztmin-ztt).GE.0._wp)
      zinbounds = MERGE(0.0_wp,zinbounds,(ztt-ztmax).GE.0._wp)
    END DO
    ! if one index was out of bounds -> print error and exit
    IF (zinbounds == 0.0_wp) CALL lookuperror(name)
    CALL fetch_ua(kidx, idx, tlucua, ua)
    IF (PRESENT(dua)) CALL fetch_ua(kidx, idx, tlucuad, dua)

  END SUBROUTINE lookup_ua_list
  !----------------------------------------------------------------------------
  SUBROUTINE lookup_uaw_list_spline(name, size, list, temp, uaw, duaw)

  ! VM: this SR is never called!

    CHARACTER(len=*),   INTENT(in)  :: name
    INTEGER,            INTENT(in)  :: size
    INTEGER,            INTENT(in)  :: list(size)
    REAL(wp),           INTENT(in)  :: temp(size)
    REAL(wp), OPTIONAL, INTENT(out) :: uaw(size), duaw(size)

    INTEGER :: idx(size)
    REAL(wp) :: zalpha(size)

    REAL(wp) :: ztt, ztshft, zinbounds, ztmax, ztmin
    INTEGER :: nl, jl

    zinbounds = 1.0_wp
    ztmin = flucupmin
    ztmax = flucupmax

    ! first compute all lookup indices and check if they are all within allowed bounds

!IBM* ASSERT(NODEPS)
    DO nl = 1, size
      jl = list(nl)
      ztshft = MERGE(1.0_wp,0.0_wp,(tmelt-temp(jl)).GE.0._wp)
      ztt = rsdeltat*temp(jl)
      zalpha(nl) = ztt - AINT(ztt)
      idx(nl) = INT(ztt-ztshft)
      zinbounds = MERGE(0.0_wp,zinbounds,(ztmin-ztt).GE.0._wp)
      zinbounds = MERGE(0.0_wp,zinbounds,(ztt-ztmax).GE.0._wp)
    END DO
    ! if one index was out of bounds -> print error and exit
    IF (zinbounds == 0.0_wp) CALL lookuperror(name, 'lookup_uaw_list_spline')
    CALL fetch_ua_spline(1, size, idx, zalpha, tlucuw, uaw, duaw)

  END SUBROUTINE lookup_uaw_list_spline
  !----------------------------------------------------------------------------
  SUBROUTINE lookup_uaw_list(name, size, list, temp, uaw, duaw)
    CHARACTER(len=*),   INTENT(in)  :: name
    INTEGER,            INTENT(in)  :: size
    INTEGER,            INTENT(in)  :: list(size)
    REAL(wp),           INTENT(in)  :: temp(size)
    REAL(wp), OPTIONAL, INTENT(out) :: uaw(size), duaw(size)

    INTEGER :: idx(size)

    REAL(wp) :: ztt, zinbounds, ztmax, ztmin
    INTEGER :: nl, jl

    zinbounds = 1.0_wp
    ztmin = fjptlucu1
    ztmax = fjptlucu2

    ! first compute all lookup indices and check if they are all within allowed bounds

!IBM* ASSERT(NODEPS)
    DO nl = 1, size
      jl = list(nl)
      ztt = ANINT(rfdeltat*temp(jl))
      idx(nl) = INT(ztt)
      zinbounds = MERGE(0.0_wp,zinbounds,(ztmin-ztt).GE.0._wp)
      zinbounds = MERGE(0.0_wp,zinbounds,(ztt-ztmax).GE.0._wp)
    END DO
    ! if one index was out of bounds -> print error and exit
    IF (zinbounds == 0.0_wp) CALL lookuperror(name, 'lookup_uaw_list')
    CALL fetch_ua(size, idx, tlucuaw, uaw)
    IF (PRESENT(duaw)) CALL fetch_ua(size, idx, tlucuaw, duaw)

  END SUBROUTINE lookup_uaw_list
  !----------------------------------------------------------------------------
  SUBROUTINE lookuperror(name, name2)
    CHARACTER(len=*), INTENT(in) :: name
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: name2

    IF (PRESENT(name2)) THEN
       WRITE(message_text, '(3a)') name, ': ', name2
    ELSE
       WRITE(message_text, '(a)') name
    END IF

    CALL finish (message_text, ' lookup table overflow (temp out of bounds)')
  END SUBROUTINE lookuperror
  !----------------------------------------------------------------------------
END MODULE mo_aes_convect_tables

