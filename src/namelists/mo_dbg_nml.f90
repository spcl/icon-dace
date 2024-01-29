! Contains the variables for debugging icon model
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

MODULE mo_dbg_nml

  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message, finish
  USE mo_impl_constants,     ONLY: max_char_length
  USE mo_io_units,           ONLY: nnml, nnml_output
  USE mo_namelist,           ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_mpi,                ONLY: my_process_is_stdio
  USE mo_nml_annotate,       ONLY: temp_defaults, temp_settings

  IMPLICIT NONE

  PUBLIC


  ! ------------------------------------------------------------------------
  ! 1.0 Namelist variables for single cell and surrounding edges/verts/neighbors
  !     debug output
  ! ------------------------------------------------------------------------

  !REAL(wp) :: s_val     = 0.0_wp  ! input  test value for salinity

  ! switches for level of debugging the icon core
  INTEGER  :: idbg_mxmn = 0       ! different levels of debug MAX/MIN output (1-5, 0: no output)
  INTEGER  :: idbg_val  = 0       ! different levels of debug output of values at lat/lon given below

  ! latitude/longitude location of single cell output for debugging
  REAL(wp) :: dbg_lat_in = 0.0_wp ! latitude of cell for debug output
  REAL(wp) :: dbg_lon_in = 0.0_wp ! longitude of cell for debug output

  ! block/index location of cell output for debugging (alternative to rlat_in/rlon_in)
  INTEGER  :: idbg_blk = 0       ! output test block
  INTEGER  :: idbg_idx = 0       ! output test index

  ! start/end vertical level for output of 3-dim variable for idbg_mxmn/val >3
  INTEGER  :: idbg_slev = -1000  ! slev=1 in sbrt dbg_mxmn
  INTEGER  :: idbg_elev =  1000  ! elev=ndim_vert in sbrt dbg_mxmn

  INTEGER, PARAMETER         :: dim_mod_tst = 20         ! array dimension of test strings
  INTEGER, PARAMETER         :: len_mod_tst = 12         ! string length of test strings
  CHARACTER(len=len_mod_tst) :: str_mod_tst(dim_mod_tst) ! namelist string of source processes to print

  NAMELIST/dbg_index_nml/  idbg_mxmn, idbg_val, dbg_lat_in, dbg_lon_in, &
    &                      idbg_blk,  idbg_idx, str_mod_tst, idbg_slev, idbg_elev

CONTAINS

  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !! Initialization of variables that control the debug output for the icon core
  !!
  !!               Initialization of single cell for debug output etc. using
  !!               namelist 'dbg_index_nml'
  !!
  SUBROUTINE read_dbg_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename
   
    INTEGER :: i_status, i
    INTEGER :: iunit
   
    CHARACTER(len=max_char_length), PARAMETER :: &
           routine = 'mo_dbg_nml/read_dbg_namelist:'

    CALL message(TRIM(routine),'read debug namelist dbg_index_nml')
   
    !------------------------------------------------------------
    ! 2.0 Read dbg_index_nml namelist
    !------------------------------------------------------------
    ! (done so far by all MPI processes)

    ! strings for marked modules to be printed out for debug purposes, exclusively
    ! default is 'all' to print out max/min and/or values in all modules

    DO i=2, dim_mod_tst
      str_mod_tst(i) = '            '
    END DO
    str_mod_tst(1) = 'all         '
   
    CALL open_nml(TRIM(filename))
    CALL position_nml ('dbg_index_nml', status=i_status)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, dbg_index_nml)  ! write defaults to temporary text file
    END IF
    SELECT CASE (i_status)
    CASE (positioned)
      READ (nnml, dbg_index_nml)                                      ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, dbg_index_nml)  ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml

    ! write the contents of the namelist to an ASCII file
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=dbg_index_nml)

  END SUBROUTINE read_dbg_namelist

END MODULE mo_dbg_nml
