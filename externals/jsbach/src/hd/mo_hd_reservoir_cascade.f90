!> Flow simulation of reservoir cascade
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2024, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
MODULE mo_hd_reservoir_cascade
!#ifndef __NO_JSBACH__
#if !defined(__NO_JSBACH__) && !defined(__NO_JSBACH_HD__)

  USE mo_kind,                   ONLY: wp

  USE mo_jsb_time,              ONLY: get_day_length

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: reservoir_cascade

  CHARACTER(len=*), PARAMETER :: modname = 'mo_hd_reservoir_cascade'

CONTAINS

  !------------------------------------------------------------------------------------------------
  SUBROUTINE reservoir_cascade(inflow, nc, steplen, retention, nres, res_content, outflow)
  !------------------------------------------------------------------------------------------------
  ! routine based on jsbach3 routine kasglob by Stefan Hagemann
  !   cascade of n equal linear reservoirs.
  !   linear reservoir approach: The outflow from a reservoir is proportional to its content.
  !------------------------------------------------------------------------------------------------

    REAL(wp),  INTENT(in)    :: inflow(:)        !< inflow to the cascade [m3/s]
    REAL(wp),  INTENT(in)    :: steplen          !< time step [s]
    INTEGER,   INTENT(in)    :: nc               !< current block size
    REAL(wp),  INTENT(in)    :: retention(:)     !< retention parameter (residence time of water in the reservoir) [day]
    REAL(wp),  INTENT(in)    :: nres(:)          !< number of reservoirs
    REAL(wp),  INTENT(inout) :: res_content(:,:) !< intermediate content of reservoir cascade [m3]
    REAL(wp),  INTENT(out)   :: outflow(:)       !< outflow at the end of the cascade [m3/s]

    ! local variables
    ! ---------------

    REAL(wp) :: &
      & day_length,        &
      & amod,              &
      & akdiv(nc),         &
      & inflow_cascade(nc)   !< inflow into each level of the cascade
    INTEGER :: nres_int(nc)  !< `nres` as integer

    INTEGER :: res, nres_max
    INTEGER :: ic

    CHARACTER(len=*), PARAMETER :: routine = modname//':reservoir_cascade'

    day_length = get_day_length()

    !$ACC DATA CREATE(nres_int, akdiv, inflow_cascade)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1) PRIVATE(amod)
    DO ic = 1, nc
      IF (nres(ic) > 0._wp) THEN ! land - without internal drainage cells
        ! unit conversions: [day] -> [seconds]
        amod  = retention(ic) * day_length
        akdiv(ic) = steplen / (amod + steplen)
      ELSE                       ! ocean and lakes
        akdiv(ic) = 0._wp
      END IF
      ! conversions: [m3/s] -> [m3]
      inflow_cascade(ic) = inflow(ic) * steplen
      !
      nres_int(ic) = NINT(nres(ic))
    END DO
    !$ACC END PARALLEL LOOP

    nres_max = SIZE(res_content(1,:))    ! maximum number of reservoirs = 'vertical' dimension of reservoir cascade

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP SEQ
    DO res = 1, nres_max
      !$ACC LOOP GANG VECTOR
      DO ic = 1, nc
        IF (res <= nres_int(ic)) THEN
          res_content(ic,res) = res_content(ic,res) + inflow_cascade(ic)
          outflow(ic) = res_content(ic,res) * akdiv(ic)
          res_content(ic,res) = res_content(ic,res) - outflow(ic)
        ELSE
          outflow(ic) = inflow_cascade(ic)
        END IF
        inflow_cascade(ic) = outflow(ic)
      END DO
      !$ACC END LOOP
    END DO
    !$ACC END LOOP
    !$ACC END PARALLEL

    ! conversions: [m3] -> [m3/s]
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO ic = 1, nc
      outflow(ic) = outflow(ic) / steplen
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC END DATA
    !$ACC WAIT(1)

  END SUBROUTINE reservoir_cascade

#endif
END MODULE mo_hd_reservoir_cascade
