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

  USE mo_time_base,              ONLY: rdaylen

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
      & amod(nc),          &
      & akdiv(nc),         &
      & inflow_cascade(nc)   !< inflow into each level of the cascade

    INTEGER :: res, nres_max

    CHARACTER(len=*), PARAMETER :: routine = modname//':reservoir_cascade'

    WHERE (nres > 0)       !< land - without internal drainage cells
      ! unit conversions: [day] -> [seconds]
      amod(:) = retention(:) * rdaylen
      akdiv(:) = 1._wp / (amod(:) + steplen) * steplen
    ELSEWHERE              !< ocean and lakes
      akdiv(:) = 0._wp
    END WHERE

    ! conversions: [m3/s] -> [m3]
    inflow_cascade(:) = inflow(:) * steplen

    nres_max = SIZE(res_content(1,:))    ! maximum number of reservoirs = 'vertical' dimension of reservoir cascade
    DO res = 1, nres_max
      WHERE (res <= NINT(nres(:)))
        res_content(:,res) = res_content(:,res) + inflow_cascade(:)
        outflow(:) = res_content(:,res) * akdiv(:)
        res_content(:,res) = res_content(:,res) - outflow(:)
      ELSEWHERE
        outflow(:) = inflow_cascade(:)
      ENDWHERE
      inflow_cascade(:) = outflow(:)
    END DO

    ! conversions: [m3] -> [m3/s]
    outflow(:) = outflow(:) / steplen

  END SUBROUTINE reservoir_cascade

#endif
END MODULE mo_hd_reservoir_cascade
