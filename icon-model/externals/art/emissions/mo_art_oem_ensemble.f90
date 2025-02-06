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

MODULE mo_art_oem_ensemble

!------------------------------------------------------------------------------
!
! Description:
!   This module contains subroutines for the generation of fluxes
!   for ensembles used in CIF. It perturbs tracers that are created by
!   OEM and are denoted as ensembles.
!
!   For the module to work, the following has to be true:
!       1. At least one tracer needs to have lensemble = .TRUE.
!       2. In the oemctrl_nml group, the parameters
!          'ens_reg_nc' (see art_oem_init_regional_map for requirements)
!          'ens_lambda_nc' (see art_oem_init_lambdas for requirements)
!       3. A working OEM - module
!       4. The number of categories and ensembles in the namelists need to
!          agree with the dimensions of the input-files. Also, the number
!          of regions in the input-files need to be consistent.
!
! This code is based on the ensemble-module written for COSMO-GHG
! by David Ochsner and adapted for ICON-ART by Michael Steiner, EMPA
!    email:  michael.steiner@empa.ch
!
!==============================================================================

  USE netcdf,                     ONLY: nf90_noerr, nf90_open, nf90_strerror, &
                                    &   nf90_inq_dimid, nf90_inq_varid,       &
                                    &   nf90_get_var, nf90_inquire_dimension, &
                                    &   nf90_nowrite, nf90_close

  USE mo_kind,                    ONLY: wp

  USE mo_exception,               ONLY: message

  USE mo_model_domain,            ONLY: t_patch, p_patch

  USE mo_parallel_config,         ONLY: nproma, idx_1d, blk_no, idx_no
  USE mo_impl_constants,          ONLY: min_rlcell,min_rlcell_int,grf_bdywidth_c
  USE mo_loopindices,             ONLY: get_indices_c
  USE mo_var_list,                ONLY: t_var_list_ptr

  !ART
  USE mo_art_data,                ONLY: p_art_data
  USE mo_art_read_xml,            ONLY: art_open_xml_file,          &
                                    &   art_close_xml_file
  USE mo_art_read_xml,            ONLY: t_xml_file
  USE mo_art_atmo_data,           ONLY: t_art_atmo
  USE mo_art_wrapper_routines,    ONLY: art_get_indices_c

  ! OEM
  USE mo_oem_config,              ONLY: ens_reg_nc,                 &
                                    &   ens_lambda_nc

  USE mo_art_oem_types,           ONLY: p_art_oem_data,             &
                                    &   t_art_oem_data,             &
                                    &   t_art_oem_config,           &
                                    &   t_art_oem_ensemble

!------------------------------------------------------------------------------

IMPLICIT NONE

  PRIVATE

  PUBLIC :: art_oem_ens_init, &
    &       art_oem_init_regional_map, &
    &       art_oem_init_lambdas

  TYPE(t_art_atmo),         POINTER :: art_atmo     !< ART atmo fields
  TYPE(t_art_oem_data),     POINTER :: oem_data     !< OEM data structure -> data
  TYPE(t_art_oem_config),   POINTER :: oem_config   !< OEM data structure -> config
  TYPE(t_art_oem_ensemble), POINTER :: oem_ensemble !< OEM data structure -> ensemble

  INTEGER, ALLOCATABLE          ::                                  &
    &  reg_map_2D(:)        ! regional map on 1D-ICON-grid DIMS:(ncell)

CONTAINS

SUBROUTINE art_oem_ens_init(jg, ierror, yerrmsg)
  !----------------------------------------------------------------------------
  ! Description:
  !   Initialize the oem_ens_init module:
  !   - check which tracers belong to the ensemble (lensemble = .TRUE.)
  !   - read in lambda values from netcdf
  !   - read in regional maps from netcdf
  !----------------------------------------------------------------------------
  INTEGER, INTENT(in)             :: jg             !< patch on which computation is performed
  INTEGER, INTENT(OUT)            :: ierror
  CHARACTER(LEN= *), INTENT(OUT)  :: yerrmsg

  INTEGER :: nt
  CHARACTER(*), PARAMETER :: routine = "oem_ens_init"

  yerrmsg = '   '
  ierror = nf90_noerr

  art_atmo => p_art_data(jg)%atmo

  oem_data => p_art_oem_data%data_fields
  oem_config => p_art_oem_data%configure
  oem_ensemble => p_art_oem_data%ensemble

  !-----------------------------------------------------------------------------
  ! Section 1: Do the setup for the module:
  !            Read regions and lambdas
  !-----------------------------------------------------------------------------
  CALL art_oem_init_regional_map(ens_reg_nc, jg, ierror, yerrmsg)
  IF (ierror /= 0) THEN
    yerrmsg = "From art_oem_init_regional_map: " // yerrmsg
    RETURN
  ENDIF

  CALL art_oem_init_lambdas(ens_lambda_nc, ierror, yerrmsg)
  IF (ierror /= 0) THEN
    yerrmsg = "From art_oem_init_lambdas: " // yerrmsg
    RETURN
  ENDIF
  ! After computing the lambda_map, we don't need the read-in data anymore
  IF (ALLOCATED(reg_map_2D)) DEALLOCATE(reg_map_2D)


END SUBROUTINE art_oem_ens_init



SUBROUTINE art_oem_init_lambdas(filepath, ierror, yerrmsg)
  !----------------------------------------------------------------------------
  ! Description:
  !   Read the lamda-values from the netcdf-file at filepath.
  !
  !   Expected structure of the NetCDF-file:
  !   
  !     Each ensemble member has a separate lambda value for each combination of region & category.
  !     Ensemble values for multiple tracers are possible.
  !
  !     The dimensions of the input-file should contain the dimensions
  !     ~(cat, reg, ens, tracer)~ and one variable:
  !
  !     1. ~lambda~ with dimensions ~(ens, reg, cat, tracer)~ and type float.
  !
  !     netcdf example_lambda_file.nc {
  !     dimensions:
  !        cat = 3 ;
  !        reg = 10 ;
  !        ens = 21 ;
  !        tracer = 2;
  !     variables:
  !        float lambda(ens, reg, cat, tracer) ;
  !               lambda:units = "unitless" ;
  !     };
  !
  !----------------------------------------------------------------------------
  CHARACTER (LEN=*), INTENT(IN)      ::                                  &
     &  filepath
  INTEGER, INTENT(OUT)               ::                                  &
     &  ierror
  CHARACTER (LEN=*), INTENT(OUT)     ::                                  &
     &  yerrmsg
  
  ! Local variables
  INTEGER                            ::                                  &
    &   ncid,         &  ! id of the openend nc file
    &   ncat,         &  ! number of categoires for ensemble
    &   nreg,         &  ! number or regions for ensemble
    &   ntra,         &  ! number of tracers for ensemble
    &   nensembles,   &  ! number of ensembles
    &   dimid,        &  ! temporary id to identify dimensions in the nc file
    &   varid,        &  ! temporary id to identify variables in the nc file
    &   i
  
  CHARACTER(*), PARAMETER :: routine="init_and_read_lambdas"
  ierror = 0
  yerrmsg = ''

  !-----------------------------------------------------------------------------
  ! Section 1: Read the lambda-values on rank 0
  !-----------------------------------------------------------------------------
  ! open the NetCDF file
  ierror = nf90_open(filepath, nf90_nowrite, ncid)
  IF (ierror /= nf90_noerr) THEN
    yerrmsg = nf90_strerror(ierror)
    RETURN
  ENDIF

  ! find size of ncat-dimension
  ierror = nf90_inq_dimid(ncid, name='cat', dimid=dimid)
  IF (ierror /= nf90_noerr) THEN
    yerrmsg = nf90_strerror(ierror)
    RETURN
  ENDIF

  ierror = nf90_inquire_dimension(ncid, dimid=dimid, len=ncat)
  IF (ierror /= nf90_noerr) THEN
    yerrmsg = nf90_strerror(ierror)
    RETURN
  ENDIF

  ! find size of nreg-dimension
  ierror = nf90_inq_dimid(ncid, name='reg', dimid=dimid)
  IF (ierror /= nf90_noerr) THEN
    yerrmsg = nf90_strerror(ierror)
    RETURN
  ENDIF

  ierror = nf90_inquire_dimension(ncid, dimid=dimid, len=nreg)
  IF (ierror /= nf90_noerr) THEN
    yerrmsg = nf90_strerror(ierror)
    RETURN
  ENDIF

  ! find size of nensembles-dimension
  ierror = nf90_inq_dimid(ncid, name='ens', dimid=dimid)
  IF (ierror /= nf90_noerr) THEN
    yerrmsg = nf90_strerror(ierror)
    RETURN
  ENDIF

  ierror = nf90_inquire_dimension(ncid, dimid=dimid, len=nensembles)
  IF (ierror /= nf90_noerr) THEN
    yerrmsg = nf90_strerror(ierror)
    RETURN
  ENDIF

  ! find size of tracer-dimension
  ierror = nf90_inq_dimid(ncid, name='tracer', dimid=dimid)
  IF (ierror /= nf90_noerr) THEN
    yerrmsg = nf90_strerror(ierror)
    RETURN
  ENDIF

  ierror = nf90_inquire_dimension(ncid, dimid=dimid, len=ntra)
  IF (ierror /= nf90_noerr) THEN
    yerrmsg = nf90_strerror(ierror)
    RETURN
  ENDIF
  ! read lambda values
  ALLOCATE(oem_data%lambda_mat(ntra, ncat, nreg, nensembles))

  ierror = nf90_inq_varid(ncid, 'lambda', varid)
  IF (ierror /= nf90_noerr) THEN
    yerrmsg = nf90_strerror(ierror)
    RETURN
  ENDIF

  ierror = nf90_get_var(ncid, varid, oem_data%lambda_mat(:,:,:,:),                  &
     &  start=(/1,1,1,1/), count=(/ntra, ncat, nreg, nensembles/))
  IF (ierror /= nf90_noerr) THEN
    yerrmsg = nf90_strerror(ierror)
    RETURN
  ENDIF

  ierror = nf90_close(ncid)
  IF (ierror /= nf90_noerr) THEN
    yerrmsg = nf90_strerror(ierror)
    RETURN
  ENDIF

END SUBROUTINE art_oem_init_lambdas


SUBROUTINE art_oem_init_regional_map(filepath, jg, ierror, yerrmsg)
  !----------------------------------------------------------------------------
  ! Description:
  !   Read the regional map from the netcdf-file at filepath
  !
  !   Expected structure of nc file:
  !
  !     The regional is defined on the 1D ICON grid.
  !     The file should contain four dimensions
  !     ~(ncell, nreg)~ and one variable:
  !
  !     2. ~REG~ with dimensions ~(reg, cell)~ and type int.
  !	       For each gridcell, have a list of integers corresponding to the
  !        different regions. Exactly one entry at each cell is expected to be
  !        1, the others are zero.
  !	       This format (opposed to one integer at each cell corresponding to a
  !        region 'index') was chosen because it's less prone to undetected
  !        errors (the number of regions is well defined because it's a dimension
  !        of the file, 'empty' regions with no associated cells can be
  !        detected easily).
  !
  !     An example of the maps-file could look like this:
  !
  !     netcdf example_reg_file.nc {
  !     dimensions:
  !        nreg  = 10 ;
  !        ncell = 84492 ;
  !     variables:
  !        int REG(reg, cell) ;
  !                REG:units = "unitless" ;
  !                REG:long_name = "region for each cell" ;
  !                REG:comment = "Exactly one entry is 1 for each
  !                               cell, the others are zero" ;
  !     }; 
  !----------------------------------------------------------------------------

  CHARACTER (LEN=*), INTENT(IN)    :: filepath
  INTEGER, INTENT(in)              :: jg      !< patch on which computation is performed
  INTEGER, INTENT(OUT)             :: ierror
  CHARACTER (LEN=*), INTENT(OUT)   :: yerrmsg

  INTEGER                          ::                                  &
    &   ncid,         &  ! id of the opened nc file
    &   nreg,         &  ! number of different regions
    &   ncell,        &  ! number of cells
    &   dimid,        &  ! temporary id to identify dimensions in the nc file
    &   cellid,       &  ! temporary id to identify dimensions in the nc file
    &   varid,        &  ! temporary id to identify variables in the nc file
    &   i, jc, jb, locidx, glbidx, i_startblk, i_endblk, is, ie, nblks_c

  CHARACTER(*), PARAMETER :: routine="art_oem_init_regional_map"

  ierror = 0
  yerrmsg = ''

  nblks_c = art_atmo%nblks

  !-----------------------------------------------------------------------------
  ! Section 1: Read the regional map
  !-----------------------------------------------------------------------------
  ierror = nf90_open(filepath, nf90_nowrite, ncid)
  IF (ierror /= nf90_noerr) THEN
    yerrmsg = nf90_strerror(ierror)
    RETURN
  ENDIF

  ierror = nf90_inq_dimid(ncid, name='cell', dimid=cellid)
  IF (ierror /= nf90_noerr) THEN
    yerrmsg = nf90_strerror(ierror)
    RETURN
  ENDIF

  ierror = nf90_inquire_dimension(ncid, cellid, len=ncell)
  IF (ierror /= nf90_noerr) THEN
    yerrmsg = nf90_strerror(ierror)
    RETURN
  ENDIF

  ALLOCATE(reg_map_2D(ncell))
  ALLOCATE(oem_data%reg_map(nproma,nblks_c))
  reg_map_2D(:) = -1
  oem_data%reg_map(:,:) = -1
    
  ierror = nf90_inq_varid(ncid, name='REG', varid=varid)
  IF (ierror /= nf90_noerr) THEN
    yerrmsg = nf90_strerror(ierror)
    RETURN
  ENDIF

  ierror = nf90_get_var(ncid, varid, reg_map_2D(:),            &
    &  start=(/1/), count=(/ncell/))
  IF (ierror /= nf90_noerr) THEN
    yerrmsg = nf90_strerror(ierror)
    RETURN
  ENDIF


  ierror = nf90_close(ncid)
  IF (ierror /= nf90_noerr) THEN
    yerrmsg = nf90_strerror(ierror)
    RETURN
  ENDIF

  !-----------------------------------------------------------------------------
  ! Section 2: Distribute the region values to the ICON internal data structure
  !-----------------------------------------------------------------------------

  i_startblk = art_atmo%i_startblk
  i_endblk = art_atmo%i_endblk
  DO jb = i_startblk, i_endblk
    ! read indices within this block:
    CALL art_get_indices_c(jg, jb, is, ie)
    DO jc = is, ie
      locidx = idx_1d(jc,jb)
      ! get global index from local index:
      glbidx = p_patch(jg)%cells%decomp_info%glb_index(locidx)
      oem_data%reg_map(jc,jb) = reg_map_2D(glbidx)
    ENDDO
  ENDDO

END SUBROUTINE art_oem_init_regional_map




FUNCTION is_lambda_map_consistent(map) RESULT(is_consistent)
  !----------------------------------------------------------------------------
  ! Description:
  !     Check if the lambdas are consistent, i.e. if all values are nonnegative.
  !----------------------------------------------------------------------------
  REAL (KIND=wp), INTENT(IN)              ::                                   &
    &  map(:,:,:,:)
  LOGICAL                                 ::                                   &
    &  is_consistent

  IF (COUNT(map < 0) > 0) THEN
    is_consistent = .FALSE.
  ELSE
    is_consistent = .TRUE.
  ENDIF

END FUNCTION is_lambda_map_consistent


!==============================================================================

END MODULE mo_art_oem_ensemble

