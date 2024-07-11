!
! mo_art_vinterp
! This module provides the vertical interpolation for external data sets
! for chemical species
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

#include "omp_definitions.inc"

MODULE mo_art_vinterp
! ICON
  USE mo_kind,                          ONLY: wp, dp
  USE mo_model_domain,                  ONLY: t_patch
  USE mo_mpi,                           ONLY: p_min, p_comm_work
  USE mo_exception,                     ONLY: message,finish,message_text
  USE mo_io_units,                      ONLY: filename_max
  USE mo_impl_constants,                ONLY: MAX_CHAR_LENGTH

  USE mo_io_config,                     ONLY: default_read_method

  USE mo_util_phys,                     ONLY: virtual_temp
  USE mo_initicon_config,               ONLY: generate_filename,  &
                                          &   ifs2icon_filename
  USE mo_master_config,                 ONLY: getModelBaseDir
  USE mo_physical_constants,            ONLY: grav, p0sl_bg
  USE mo_nh_vert_interp,                ONLY: prepare_extrap
  USE mo_nh_vert_interp_ipz,            ONLY: z_at_plevels
  USE mo_read_interface,                ONLY: t_stream_id, nf, openInputFile,  &
                                          &   closeFile,  read_2d_1time,       &
                                          &   read_2d_1lev_1time,              &
                                          &   read_3d_1time, on_cells

  USE mo_initicon_types,                ONLY: geop_ml_var
! ART
  USE mo_art_data,                      ONLY: p_art_data
  USE mo_art_atmo_data,                 ONLY: t_art_atmo
  USE mo_art_chem_init_types,           ONLY: t_art_chem_init_coord, &
                                          &   t_chem_init_state
  USE mo_art_chem_init_coord,           ONLY: alloc_vct_chem_init,           &
                                          &   init_vct_chem_init,            &
                                          &   half_level_pressure_chem_init, &
                                          &   auxhyb_chem_init,              &
                                          &   full_level_pressure_chem_init, & 
                                          &   geopot_chem_init
  USE mo_art_chem_init_utils,           ONLY: allocate_chem_init_atm,  &
                                          &   allocate_chem_init_chem, &
                                          &   deallocate_chem_init_atm
  IMPLICIT NONE 
  
  PRIVATE


  PUBLIC  :: art_prepare_vinterp
  PUBLIC  :: art_prepare_vinterp_pres
  
  INCLUDE 'netcdf.inc'
  
CONTAINS


SUBROUTINE art_prepare_vinterp(chem_init,p_patch, &
                         &     coord_path, chem_init_path)
!<
! SUBROUTINE art_prepare_vinterp_pres
! This subroutine computes the geopotential height of the levels
! of a dataset on hybrid levels
! Part of Module: mo_art_vinterp
! Author: Jennifer Schroeter, KIT
! Initial Release: 2015-11-15
! Modifications:
! 2017-09-06: Michael Weimer, KIT
! - generalised algorithm to interpolate external datasets (not only
!   initialisation)
!>

  TYPE(t_chem_init_state), INTENT(inout), TARGET  :: &
    &  chem_init              !< container with meta for vertical interpolation

  TYPE(t_patch), TARGET,INTENT(IN) ::  & 
    &  p_patch               !< patch on which computation is performed

  CHARACTER(LEN=*), INTENT(IN) ::  &
    &  coord_path             !< file name of coord file  
                              !  (info about vertical coordinate)
  CHARACTER(LEN=*), INTENT(IN) ::  &
    &  chem_init_path         !< file name of the file with which 
                              !  the tracer is prescribed

  ! local variables
  INTEGER :: ncid, dimid, varid,  varid_geop    !< netCDF variables

  CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
    &  routine = 'art_prepare_vininterp'

  INTEGER ::      &
    &  no_levels, &  !< number of levels in the external dataset
    &  jg,        &  !< patch id
    &  jb,        &  !< loop index
    &  nlen          !< filled values in the first horizontal dimension

  REAL(wp), ALLOCATABLE ::  &
    &  pres_ic_chem_init(:,:)  !< half level pressures in the external data (Pa)
  REAL(wp), ALLOCATABLE ::     &
    &  delp_chem_init(:,:),    & !< arrays for vertical interpolation
    &  rdelp_chem_init(:,:),   & 
    &  rdlnpr_chem_init(:,:),  &
    &  rdalpha_chem_init(:,:), &
    &  geop_mc_chem_init(:,:)    !< geopotential at full levels
  REAL(wp), ALLOCATABLE ::   &
    &  lnp_ic_chem_init(:,:)     !< logarithmic pressure ?
  REAL(wp), ALLOCATABLE ::  &
    &  geop_ic_chem_init(:,:)    !< geopotential at interfaces

  CHARACTER(LEN=filename_max) ::  &
    &  ifs2icon_file     !< filename of the IFS initialisation file
  CHARACTER(LEN=10) ::  &
    &  geop_ml_var_chem_init      !< name of the surface geopotential variable 
                                  !  in external file

  REAL(wp) ::  &
    &  min_tv, local_min_tv  !< minimum virtual temperature 
                             !  globally and locally for this process

  LOGICAL ::             &
    &  l_exist,          & !< flag if IFS file exists
    &  lgeop_sfc_in_file   !< flag if surface geopotential exists 
                           !  in external dataset

  TYPE(t_art_atmo), POINTER :: &
    &  art_atmo            !< ART atmo fields
  TYPE(t_art_chem_init_coord), POINTER :: &
    &  chem_coord          !< for abbreviation of chem_init%chem_init_coord only
  TYPE(t_stream_id) ::  &
    &  stream_id
 
  INTEGER ::  &
    &  geop_sfc_var_ndims  !< number of dimensions in geopotential variable
  CHARACTER(LEN=10) ::  &
    &  geop_sfc_var !< surface-level surface geopotential
  REAL(wp), ALLOCATABLE :: vct_chem_init_a(:) ! param. A of the vertical coordinte
  REAL(wp), ALLOCATABLE :: vct_chem_init_b(:) ! param. B of the vertical coordinate

  jg = p_patch%id

  art_atmo => p_art_data(jg)%atmo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Check for existence of necessary files for vertical interpolation
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Check if cart_cheminit_file exists
  INQUIRE(file = TRIM(chem_init_path), EXIST = l_exist)
  IF (.NOT. l_exist) THEN
    CALL finish('mo_art_vinterp:art_prepare_vinterp',                    &
          &     'cart_cheminit_file '//TRIM(chem_init_path)              &
          &   //' not found.')
  END IF

  ! Check if cart_cheminit_coord exists
  INQUIRE(file = TRIM(coord_path), EXIST = l_exist)
  IF (.NOT. l_exist) THEN
    CALL finish('mo_art_vinterp:art_prepare_vinterp',                    &
          &     'cart_cheminit_coord '//TRIM(coord_path)                 &
          &   //' not found.')
  END IF
  
  ! Check if IFS file exists
  ifs2icon_file = generate_filename(ifs2icon_filename, getModelBaseDir(), &
                      &         art_atmo%nroot,art_atmo%nbisect, jg)
  INQUIRE(file = TRIM(ifs2icon_file), EXIST = l_exist)

  IF (.NOT. l_exist) THEN
    CALL finish('mo_art_vinterp:art_prepare_vinterp',                   &
          &     'IFS file '//TRIM(ifs2icon_file)                        &
          &   //' not found. If this is a nested run please consider '  &
          &   //'setting grid_nml:start_time greater than zero for the '&
          &   //'nest (e.g. twice of the time step).')
  END IF


  chem_coord => chem_init%chem_init_coord
 
  CALL nf(nf_open(TRIM(coord_path), NF_NOWRITE, ncid), routine)


  CALL nf(nf_inq_dimid(ncid, 'lev', dimid), routine)
  CALL nf(nf_inq_dimlen(ncid, dimid, no_levels), routine)

  chem_init%chem_init_in%nlev_chem_init = no_levels

  ALLOCATE(chem_coord%vct_chem_init(2*(no_levels+1)))

  ALLOCATE(vct_chem_init_a(no_levels+1))
  ALLOCATE(vct_chem_init_b(no_levels+1))
  ALLOCATE(pres_ic_chem_init(art_atmo%nproma, no_levels+1))
  ALLOCATE(lnp_ic_chem_init(art_atmo%nproma, no_levels+1))
  ALLOCATE(delp_chem_init(art_atmo%nproma, no_levels))
  ALLOCATE(rdelp_chem_init(art_atmo%nproma, no_levels))
  ALLOCATE(rdlnpr_chem_init(art_atmo%nproma, no_levels))
  ALLOCATE(rdalpha_chem_init(art_atmo%nproma, no_levels))
  ALLOCATE(geop_mc_chem_init(art_atmo%nproma, no_levels))
  ALLOCATE(geop_ic_chem_init(art_atmo%nproma, no_levels+1))

  IF (nf_inq_varid(ncid, 'hyai', varid) == nf_noerr) THEN
    CALL nf(nf_get_var_double(ncid, varid, vct_chem_init_a), routine)  
  ELSE
    CALL finish(TRIM(routine),' hyai is missing')
  ENDIF

  IF (nf_inq_varid(ncid, 'hybi', varid) == nf_noerr) THEN
    CALL nf(nf_get_var_double(ncid, varid, vct_chem_init_b), routine)
  ELSE
    CALL finish(TRIM(routine),' hybi is missing')
  ENDIF

  CALL nf(nf_close(ncid), routine)

  SELECT CASE (TRIM(ADJUSTL(chem_init%chem_init_in%model_name)))
    CASE ('MOZART')
      CALL message(routine, "THIS IS MOZART INIT")
      chem_coord%vct_chem_init(1:no_levels+1)  =   p0sl_bg * vct_chem_init_a(:)
    CASE ('EMAC', 'ERA')
      chem_coord%vct_chem_init(1:no_levels+1)  =   vct_chem_init_a(:)
    CASE DEFAULT
      CALL message('','=============================================')
      CALL message('','WARNING: unknown model name for vertical interpolation. '  &
          &         //'Assuming EMAC configuration: press = hyai + hybi * PS '    &
          &         //'with hyai given in Pa.')
      CALL message('','=============================================')

      chem_coord%vct_chem_init(1:no_levels+1)             =   vct_chem_init_a(:)
  END SELECT

  ! vct_chem_init consists of the hyai parameter in the first half and hybi in
  ! the second half of the array
  chem_coord%vct_chem_init(no_levels+2:2*(no_levels+1)) = vct_chem_init_b(:)

  CALL alloc_vct_chem_init(chem_coord,no_levels)
  CALL init_vct_chem_init(chem_coord,no_levels)
  CALL allocate_chem_init_atm(chem_init, jg ,no_levels)
  CALL allocate_chem_init_chem(chem_init,jg,no_levels)


  ! try to get surface geopotential from chem_init file if it is available 
  ! (else read it from IFS initialisation)
  ! Check additionally for necessary variables T, Q, and PS


  CALL nf(nf_open(TRIM(chem_init_path), NF_NOWRITE, ncid), routine)

  IF (nf_inq_varid(ncid, 'T', varid) /= nf_noerr) THEN
    CALL finish(TRIM(routine),'T is missing')
  ENDIF

  IF (nf_inq_varid(ncid, 'Q', varid) /= nf_noerr) THEN
    CALL finish(TRIM(routine),'Q is missing')
  ENDIF

  IF (nf_inq_varid(ncid, 'PS', varid) /= nf_noerr) THEN
    CALL finish(TRIM(routine),'PS is missing')
  ENDIF

  CALL nf(nf_close(ncid), routine)


  
  CALL nf(nf_open(TRIM(ifs2icon_file), NF_NOWRITE, ncid), routine)

  IF (nf_inq_varid(ncid, 'GEOSP', varid_geop) == nf_noerr) THEN
    geop_ml_var_chem_init = 'GEOSP'
  ELSE IF (nf_inq_varid(ncid, 'GEOP_ML', varid_geop) == nf_noerr) THEN
    geop_ml_var_chem_init= 'GEOP_ML'
  ELSE
    CALL finish(TRIM(routine),'Could not find model-level sfc geopotential')
  ENDIF

  IF (.NOT. (nf_inq_varid(ncid, TRIM(geop_ml_var_chem_init), varid_geop) == nf_noerr)) THEN
    CALL finish(TRIM(routine),'surface geopotential var '   &
                       &     //TRIM(geop_ml_var_chem_init)//' is missing')
  ENDIF

  CALL openInputFile(stream_id, ifs2icon_file, p_patch, &
    &                default_read_method)

  lgeop_sfc_in_file = (nf_inq_varid(ncid, 'GEOP_SFC', varid) == nf_noerr)

  IF (lgeop_sfc_in_file) THEN
    geop_sfc_var = 'GEOP_SFC'
  ELSE

    WRITE (message_text,'(a,a)')                            &
      &  'surface-level surface geopotential is missing. ', &
      &  'use model-level surface geopotential, instead.'
    CALL message(TRIM(routine),TRIM(message_text))

    ! use model level geopotential instead
    geop_sfc_var = geop_ml_var
  ENDIF


  IF (nf_inq_varid(ncid, TRIM(geop_sfc_var), varid) == nf_noerr) THEN
    CALL nf(nf_inq_varndims(ncid,varid,geop_sfc_var_ndims), routine)
  ELSE
    CALL finish(TRIM(ROUTINE),'surface geopotential variable '//TRIM(geop_sfc_var)//' is missing')
  ENDIF

  CALL nf(nf_close(ncid), routine)
  
  IF (lgeop_sfc_in_file) THEN

    IF (geop_sfc_var_ndims == 3) THEN
      CALL read_2d_1lev_1time(stream_id, on_cells, TRIM(geop_sfc_var),   &
        &                     fill_array= chem_init%chem_init_in%phi_sfc)
    ELSE IF (geop_sfc_var_ndims == 2) THEN
      message_text = 'reading in GEOSP in 2d. '
      CALL message(TRIM(routine),TRIM(message_text))

      CALL read_2d_1time(stream_id, on_cells, TRIM(geop_sfc_var), &
        &                     fill_array= chem_init%chem_init_in%phi_sfc)
    ELSE
      CALL finish(TRIM(routine),"geop_sfc_var: Dimension mismatch")
    ENDIF
    
  END IF

  CALL closeFile(stream_id)


  CALL openInputFile( stream_id, TRIM(chem_init_path), p_patch, &
    &                 default_read_method)
  CALL read_3d_1time(stream_id, on_cells,'T', fill_array= chem_init%chem_init_in%temp)
  CALL read_3d_1time(stream_id, on_cells,'Q', fill_array= chem_init%chem_init_in%q)
  CALL read_2d_1time(stream_id, on_cells,'PS', fill_array= chem_init%chem_init_in%ps)

  CALL nf(nf_open(TRIM(chem_init_path), NF_NOWRITE, ncid), routine)
  

  ! if virtual temperature is in the dataset itself, take it, else calculate it
  ! internally with temperature and water vapour of the external file
  IF (nf_inq_varid(ncid, 'TV', varid) == nf_noerr) THEN
    CALL read_3d_1time(stream_id, on_cells,'TV', fill_array= chem_init%chem_init_in%temp_v)

    ! compute global minimum of virtual temperature and test
    local_min_tv = MINVAL(chem_init%chem_init_in%temp_v(:,:,1:art_atmo%nblks-1))
    min_tv = p_min(REAL(local_min_tv, dp),comm = p_comm_work)
    local_min_tv = min_tv
    
    IF (min_tv  < 0._wp) THEN
      CALL message(routine,'TV is negative in '//TRIM(chem_init_path)  &
        &                //'. Trying to calculate it by variables T and Q...')

      CALL virtual_temp(p_patch,chem_init%chem_init_in%temp,        &
               &        chem_init%chem_init_in%q,                   &
               &        temp_v=chem_init%chem_init_in%temp_v)
    END IF
          
  ELSE
    CALL virtual_temp(p_patch,chem_init%chem_init_in%temp,        &
             &        chem_init%chem_init_in%q,                   &
             &        temp_v=chem_init%chem_init_in%temp_v)
  END IF


  CALL nf(nf_close(ncid), routine)
  CALL closeFile(stream_id)


  DO jb = 1,art_atmo%nblks
    IF (jb /= art_atmo%nblks) THEN
       nlen = art_atmo%nproma
    ELSE
       nlen = art_atmo%npromz
    ENDIF

    IF (MAXVAL(chem_init%chem_init_in%ps(1:nlen,jb)) <= 100._wp) THEN
      chem_init%chem_init_in%ps(1:nlen,jb)    &
          &       = EXP(chem_init%chem_init_in%ps(1:nlen,jb))
    ENDIF

    CALL half_level_pressure_chem_init(chem_coord,chem_init%chem_init_in%ps(:,jb),   &
                           &           art_atmo%nproma, nlen, no_levels, pres_ic_chem_init)

    CALL full_level_pressure_chem_init(chem_coord,pres_ic_chem_init,art_atmo%nproma,   &
                           &           nlen, no_levels,                                &
                           &           chem_init%chem_init_in%pres(:,:,jb))

    CALL auxhyb_chem_init(chem_coord,pres_ic_chem_init, art_atmo%nproma, nlen,    &
                 &        no_levels,delp_chem_init, rdelp_chem_init,              &
                 &        lnp_ic_chem_init, rdlnpr_chem_init, rdalpha_chem_init) 

    CALL geopot_chem_init(chem_coord,chem_init%chem_init_in%temp_v(:,:,jb),              &
                 &        rdlnpr_chem_init, rdalpha_chem_init,                           &
                 &        chem_init%chem_init_in%phi_sfc(:,jb),art_atmo%nproma,1, nlen,  &
                 &        no_levels, geop_mc_chem_init, geop_ic_chem_init )
  
    chem_init%chem_init_in%z3d(1:nlen,1:no_levels,jb)  &
             &      = geop_mc_chem_init(1:nlen,1:no_levels) / grav

  ENDDO



END SUBROUTINE art_prepare_vinterp
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_prepare_vinterp_pres(jg,chem_init,file_namepath)
!<
! SUBROUTINE art_prepare_vinterp_pres
! This subroutine computes the height of the pressure levels
! of a dataset on pressure levels (basically by using z_at_plevels)
! ATTENTION: currently this routines assumes pressure values in units of hPa,
!            not SI units
! Part of Module: mo_art_vinterp
! Author: Michael Weimer, KIT
! Initial Release: 2017-08-11
! Modifications:
!>
  IMPLICIT NONE
  INTEGER, INTENT(in) :: &
    &  jg
  TYPE(t_chem_init_state), INTENT(inout) :: &
    &  chem_init         !< structure with meta data of the dataset
  CHARACTER(LEN=*), INTENT(in) :: &
    &  file_namepath     !< file name and path of the dataset
  !local variables
  INTEGER :: &
    &  no_lev                  !< number of levels in the dataset
  INTEGER :: jc, jk,jb         !< loop indices
  TYPE(t_art_atmo), POINTER :: &
    &  art_atmo               !< ART atmo fields

  art_atmo => p_art_data(jg)%atmo
  
  ! determine number of model levels in the external data
  IF (.NOT. chem_init%chem_init_in%linitialized) THEN
    CALL art_get_vertical_dimension(TRIM(file_namepath),no_lev)

    ALLOCATE(chem_init%chem_init_in%z3d(art_atmo%nproma,no_lev,art_atmo%nblks))
    chem_init%chem_init_in%nlev_chem_init = no_lev
    CALL allocate_chem_init_chem(chem_init,jg,no_lev)

    ALLOCATE(chem_init%chem_init_in%pres_pl(no_lev))
    ALLOCATE(chem_init%chem_init_in%pres_pl3d(art_atmo%nproma,no_lev,art_atmo%nblks))
    ALLOCATE(chem_init%chem_init_in%kpbl1(art_atmo%nproma,art_atmo%nblks))
    ALLOCATE(chem_init%chem_init_in%kpbl2(art_atmo%nproma,art_atmo%nblks))
    ALLOCATE(chem_init%chem_init_in%wfacpbl1(art_atmo%nproma,art_atmo%nblks))
    ALLOCATE(chem_init%chem_init_in%wfacpbl2(art_atmo%nproma,art_atmo%nblks))

    chem_init%chem_init_in%linitialized = .TRUE.
  ELSE 
    no_lev = chem_init%chem_init_in%nlev_chem_init
  END IF


  ! read pressure coordinate
  CALL art_read_lev_variable(TRIM(file_namepath),chem_init%chem_init_in%pres_pl)
  chem_init%chem_init_in%pres_pl = chem_init%chem_init_in%pres_pl * 100._wp  ! because of hPa

  ! make it a three dimensional array
  DO jb = 1, art_atmo%nblks
    DO jk = 1,no_lev
      DO jc = 1, art_atmo%nproma
        chem_init%chem_init_in%pres_pl3d(jc,jk,jb) = chem_init%chem_init_in%pres_pl(jk)
      END DO
    END DO
  END DO

  ! get factors for vertical interpolation
  CALL prepare_extrap(art_atmo%z_mc,                    & !in
              &       art_atmo%nblks,                   & !in
              &       art_atmo%npromz,                  & !in
              &       art_atmo%nlev,                    & !in
              &       chem_init%chem_init_in%kpbl1,     & !out
              &       chem_init%chem_init_in%wfacpbl1,  & !out
              &       chem_init%chem_init_in%kpbl2,     & !out
              &       chem_init%chem_init_in%wfacpbl2 )   !out

  
  CALL z_at_plevels(art_atmo%pres, art_atmo%tempv,      & !in
            &       art_atmo%z_mc,                      & !in
            &       chem_init%chem_init_in%pres_pl3d,   & !in
            &       chem_init%chem_init_in%z3d,         & !out
            &       art_atmo%nblks,art_atmo%npromz,     & !in
            &       art_atmo%nlev, no_lev,              & !in
            &       chem_init%chem_init_in%kpbl1,       & !in
            &       chem_init%chem_init_in%wfacpbl1,    & !in
            &       chem_init%chem_init_in%kpbl2,       & !in
            &       chem_init%chem_init_in%wfacpbl2)      !in

END SUBROUTINE art_prepare_vinterp_pres
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_get_vertical_dimension(pathfile,dim_levs)
!<
! SUBROUTINE art_get_vertical_dimension
! This subroutine reads the number of levels in the dataset
! Part of Module: mo_art_vinterp
! Author: Michael Weimer, KIT
! Initial Release: about 2017-05-12
! Modifications:
!>
  IMPLICIT NONE
  CHARACTER(LEN = *), INTENT(in)  :: &
    &  pathfile
  INTEGER, INTENT(out)  ::  &
    &  dim_levs

  ! local variables
  CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
    routine = 'mo_art_read_emissions:art_get_vertical_dimension'
  INTEGER :: &
    & ncid, varid, dimids_lev(1)

  CALL nf(nf_open(TRIM(pathfile), NF_NOWRITE, ncid), routine)

  IF (nf_inq_varid(ncid, 'lev', varid) == nf_noerr) THEN
    CALL nf(nf_inq_vardimid(ncid,varid,dimids_lev), routine)
    CALL nf(nf_inq_dimlen(ncid,dimids_lev(1),dim_levs), routine)
  ELSE
    CALL finish(TRIM(routine),'Could not find "lev" variable in '  &
                     &      //TRIM(pathfile))
  ENDIF
  CALL nf(nf_close(ncid), routine)
END SUBROUTINE art_get_vertical_dimension
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_read_lev_variable(pathfile,levs)
!<
! SUBROUTINE art_read_lev_variable
! This subroutine opens, reads and closes a committed dataset and reads 
! the variable "lev"
! Part of Module: mo_art_vinterp
! Author: Michael Weimer, KIT
! Initial Release: about 2016-11-17
! Modifications:
!>
  IMPLICIT NONE
  CHARACTER(LEN = *), INTENT(in)  :: &
    &  pathfile
  REAL(wp), INTENT(out)  ::  &
    &  levs(:)
  ! local variables
  CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
    routine = 'mo_art_read_emissions:art_read_lev_variable'
  INTEGER :: &
    & ncid, varid

  CALL nf(nf_open(TRIM(pathfile), NF_NOWRITE, ncid), routine)

  IF (nf_inq_varid(ncid, 'lev', varid) == nf_noerr) THEN
    CALL nf(nf_get_var_double(ncid, varid, levs), routine)
  ELSE
    CALL finish(TRIM(routine),'Could not find "lev" variable in '  &
                     &      //TRIM(pathfile))
  ENDIF
  CALL nf(nf_close(ncid), routine)
END SUBROUTINE art_read_lev_variable
!!
!!-------------------------------------------------------------------------
!!
END MODULE


