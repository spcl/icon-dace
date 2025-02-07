! mo_jsbach_comm_to_echam5mods.f90 - Domain characteristics for JSBACH usage
!
! Copyright (C) 2014, MPI-M
! SPDX-License-Identifier: BSD-3-Clause
! See ./LICENSES/ for license information
!_________________________________________

MODULE mo_jsbach_comm_to_echam5mods
  ! 
  ! Description: 
  !   <Say what this module is for> 
  ! 
  ! Current Code Owner: <Name of person responsible for this code> 
  ! 
  ! History: 
  !  
  ! Version   Date     Comment 
  ! -------   ----     ------- 
  ! <version> <date>   Original code. <Your name> 
  ! 
  ! Code Description: 
  !   Language:           Fortran 90. 
  !   Software Standards: "European Standards for Writing and  
  !     Documenting Exchangeable Fortran 90 Code". 
  ! 
  !  This routine originates (year 2014) from MPI-ESM, the Earth System Model of the 
  !  Max Planck Institute for Meteorology (Mauritsen et al. 2019). 
  !  Reference: Mauritsen, T., et al. (2019) Developments in the MPI-M Earth System Model 
  !  version 1.2 (MPI-ESM1.2) and its response to increasing CO2. J. Adv. Model. Earth Syst., 11, 
  !  doi: 10.1029/2018MS001400.

  USE mo_kind, ONLY: dp

  IMPLICIT NONE 

  PUBLIC

  INTEGER           :: nlon, nlat, nland
  INTEGER,  POINTER :: kpoints(:)
  REAL(dp), POINTER :: lon(:), lat(:)
  LOGICAL,  POINTER :: mask(:,:)

  INTEGER           :: domain_nlon, domain_nlat, domain_nland
  LOGICAL,  POINTER :: domain_mask(:,:)

END MODULE mo_jsbach_comm_to_echam5mods

!- End of module header
