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

MODULE mo_reader_abstract

  USE mo_kind,         ONLY: dp
  USE mo_model_domain, ONLY: t_patch
  USE mtime,           ONLY: julianday
  
  IMPLICIT NONE
  
  PUBLIC t_abstract_reader

  TYPE, ABSTRACT :: t_abstract_reader
  CONTAINS
    PROCEDURE(abstract_init),            DEFERRED :: init
    PROCEDURE(abstract_get_one_timelev), DEFERRED :: get_one_timelev
    PROCEDURE(abstract_deinit),          DEFERRED :: deinit
    PROCEDURE(abstract_get_times),       DEFERRED :: get_times
    ! we need those, since we only assume the read data to be blocked
    ! we make no assumptions for the data to be placed on cells,
    ! edges, lonlat grids, ...
    PROCEDURE(abstract_get_nblks),       DEFERRED :: get_nblks
    PROCEDURE(abstract_get_npromz),      DEFERRED :: get_npromz
  END TYPE t_abstract_reader

  ABSTRACT INTERFACE
    SUBROUTINE abstract_init (this, p_patch, filename)
      IMPORT :: t_abstract_reader, t_patch
      CLASS(t_abstract_reader), INTENT(inout) :: this
      TYPE(t_patch),    TARGET, INTENT(in   ) :: p_patch
      CHARACTER(len=*),         INTENT(in   ) :: filename
    END SUBROUTINE abstract_init

    SUBROUTINE abstract_get_one_timelev (this, timelevel, varname, dat)
      IMPORT :: t_abstract_reader, dp
      CLASS(t_abstract_reader), INTENT(inout) :: this
      INTEGER,                  INTENT(in   ) :: timelevel
      CHARACTER(len=*),         INTENT(in   ) :: varname
      REAL(dp), ALLOCATABLE,    INTENT(inout) :: dat(:,:,:,:)
    END SUBROUTINE abstract_get_one_timelev

    SUBROUTINE abstract_get_times (this, times)
      IMPORT :: t_abstract_reader, julianday
      CLASS(t_abstract_reader),     INTENT(inout) :: this
      TYPE(julianday), ALLOCATABLE, INTENT(  out) :: times(:)
    END SUBROUTINE abstract_get_times

    SUBROUTINE abstract_deinit (this)
      IMPORT :: t_abstract_reader
      CLASS(t_abstract_reader), INTENT(inout) :: this
    END SUBROUTINE abstract_deinit

    FUNCTION abstract_get_nblks (this) RESULT(nblks)
      IMPORT :: t_abstract_reader
      CLASS(t_abstract_reader), INTENT(in   ) :: this
      INTEGER                                 :: nblks
    END FUNCTION

    FUNCTION abstract_get_npromz (this) RESULT(npromz)
      IMPORT :: t_abstract_reader
      CLASS(t_abstract_reader), INTENT(in   ) :: this
      INTEGER                                 :: npromz
    END FUNCTION
  END INTERFACE

  CHARACTER(len=*), PARAMETER :: modname = 'mo_abstract_reader'

END MODULE mo_reader_abstract
