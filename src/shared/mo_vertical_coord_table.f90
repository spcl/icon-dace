! module *mo_vertical_coord_table* - *loop indices and surface-pressure independent
! variables associated with the vertical finite-difference scheme.
!
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

MODULE mo_vertical_coord_table

  !-------------------------------------------------------------------------
  !
  !    ProTeX FORTRAN source: Style 2
  !    modified for ICON project, DWD/MPI-M 2006
  !
  !-------------------------------------------------------------------------
  !
  !
  !

  ! USE mo_parameters

  USE mo_kind,               ONLY: wp
  USE mo_io_units,           ONLY: filename_max, find_next_free_unit
  USE mo_exception,          ONLY: message_text, message, finish
  USE mo_impl_constants,     ONLY: SUCCESS, max_char_length

  IMPLICIT NONE

  PRIVATE


  PUBLIC :: allocate_vct_atmo
  PUBLIC :: read_vct
  PUBLIC :: vct_a, vct_b, vct


  REAL(wp), ALLOCATABLE :: vct_a(:) ! param. A of the vertical coordinte
  REAL(wp), ALLOCATABLE :: vct_b(:) ! param. B of the vertical coordinate
  REAL(wp), ALLOCATABLE :: vct  (:) ! param. A and B of the vertical coordinate


  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_vertical_coord_table'


CONTAINS

  !EOC
  !-------------------------------------------------------------------------

  !>
  !! Read the A and B parameters of the hybrid vertical grid.
  !!
  !! Read the A and B parameters of the hybrid vertical grid,
  !! which define the half level heights: zh=A+B*topo [m]
  !!
  SUBROUTINE read_vct (klev, vct_file, vct_a, vct_b)

    INTEGER,           INTENT(IN   ) :: klev
    CHARACTER(LEN=*),  INTENT(IN   ) :: vct_file
    REAL(wp),          INTENT(INOUT) :: vct_a(:), vct_b(:)  

    ! Local variables
    CHARACTER(len=max_char_length),PARAMETER :: routine  = &
         &   'mo_vertical_coord_table:read_vct'
    CHARACTER(len=filename_max)              :: line

    INTEGER :: ist, iunit
    INTEGER :: ik, jk

    !-------------------------------------------------------------------------
    !BOC


    ! Open file
    !
    ! use hybrid sigma height tables
    iunit = find_next_free_unit(10,20)
    OPEN (unit=iunit,file=TRIM(vct_file),access='SEQUENTIAL', &
      &  form='FORMATTED', action='READ', status='OLD', IOSTAT=ist)

    IF(ist/=success)THEN
      CALL finish (TRIM(routine), 'open vertical coordinate table file failed')
    ENDIF

    ! Skip header line
    READ (iunit,*,IOSTAT=ist) line
    IF(ist/=success)THEN
      CALL finish (TRIM(routine), 'reading header line failed')
    ENDIF

    ! Read A and B
    DO jk=1,klev+1
       READ (iunit,*,IOSTAT=ist) ik, vct_a(jk), vct_b(jk)
       IF(ist/=success)THEN
          CALL finish (TRIM(routine), 'reading vct_a and vct_b failed')
       ENDIF
    END DO
    !$ACC UPDATE DEVICE(vct_a) ASYNC(1)

    CALL message(TRIM(routine), 'vertical coordinate table file successfully read')

    CLOSE(iunit)

  END SUBROUTINE read_vct

  !EOC
  !-------------------------------------------------------------------------


  !----------------------------------------------------------------------------------------------------
  !> Utility function: Allocation of vertical coordinate tables.
  !!
  SUBROUTINE allocate_vct_atmo(nlevp1)
    INTEGER,                INTENT(IN)    :: nlevp1
    ! local variables
    CHARACTER(*), PARAMETER   :: routine = modname//"::allocate_vct_atmo"
    INTEGER :: error_status

    ! Allocate input for init routines
    ALLOCATE(vct_a(nlevp1), vct_b(nlevp1), STAT=error_status)
    !$ACC ENTER DATA CREATE(vct_a)
    IF (error_status/=SUCCESS) CALL finish (TRIM(routine), 'allocation of vct_a/vct_b failed')
    
    ! Allocate input for derived variables of input
    ALLOCATE(vct(nlevp1*2), STAT=error_status)
    IF (error_status/=SUCCESS) CALL finish (TRIM(routine), 'allocation of vct failed')
    
  END SUBROUTINE allocate_vct_atmo

END MODULE mo_vertical_coord_table

