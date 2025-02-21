! Namelist reading for synthetic radar data on the model grid
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

MODULE mo_synradar_nml


  IMPLICIT NONE
  PUBLIC :: read_synradar_namelist

  ! module name
  CHARACTER(*), PARAMETER :: modname = "mo_synradar_nml"
  
CONTAINS
  !! Read Namelist for I/O.
  !!
  !! This subroutine
  !! - reads the Namelist for I/O
  !! - sets default values
  !! - potentially overwrites the defaults by values used in a
  !!   previous integration (if this is a resumed run)
  !! - reads the user's (new) specifications
  !! - stores the Namelist for restart
  !! - fills the configuration state (partly)
  !!
  SUBROUTINE read_synradar_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN)   :: filename

    
  END SUBROUTINE read_synradar_namelist

END MODULE mo_synradar_nml
