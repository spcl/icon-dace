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

MODULE mo_extpar_config

  USE mo_kind,               ONLY: wp
  USE mo_impl_constants,     ONLY: max_dom, MAX_CHAR_LENGTH
  USE mo_io_units,           ONLY: filename_max
  USE mo_util_string,        ONLY: t_keyword_list, &
    &                              associate_keyword, with_keywords, &
    &                              int2string
  USE mo_exception,          ONLY: finish


  IMPLICIT NONE

  PRIVATE

  ! variables
  PUBLIC :: itopo
  PUBLIC :: fac_smooth_topo
  PUBLIC :: n_iter_smooth_topo
  PUBLIC :: hgtdiff_max_smooth_topo
  PUBLIC :: itype_lwemiss
  PUBLIC :: read_nc_via_cdi
  PUBLIC :: heightdiff_threshold
  PUBLIC :: lrevert_sea_height
  PUBLIC :: pp_sso
  PUBLIC :: itype_vegetation_cycle
  PUBLIC :: extpar_filename
  PUBLIC :: extpar_varnames_map_file
  PUBLIC :: i_lctype
  PUBLIC :: nclass_lu
  PUBLIC :: nhori          ! number of sectors for horizon
  PUBLIC :: nmonths_ext

  ! subroutines/functions
  PUBLIC :: generate_filename
  PUBLIC :: generate_td_filename

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_extpar_config'


  !>
  !!----------------------------------------------------------------------------
  !! Derived type containing control variables specific to the nonhydrostatic 
  !! atm model

  !------------------------------------------------------------------------

  ! namelist variables

  INTEGER  :: itopo       ! 0: topography specified by analytical functions,
                          ! 1: topography read from netcdf files

  REAL(wp) :: fac_smooth_topo
  INTEGER  :: n_iter_smooth_topo(max_dom)
  REAL(wp) :: hgtdiff_max_smooth_topo(max_dom)
  INTEGER  :: itype_lwemiss     ! switch to select longwave emissivity data
  LOGICAL  :: read_nc_via_cdi ! read netcdf input via cdi library (alternative: parallel netcdf)
  REAL(wp) :: heightdiff_threshold(max_dom)
  LOGICAL  :: lrevert_sea_height  ! if true: bring sea points back to original height
  INTEGER  :: pp_sso      ! if >0: postprocess SSO over glaciers to reduce contribution of mean slope
  INTEGER  :: itype_vegetation_cycle
  INTEGER  :: nhori           ! number of sectors for horizon
  !$ACC DECLARE COPYIN(nhori)

  ! ExtPar input filename, may contain keywords, by default
  ! extpar_filename = "<path>extpar_<gridfile>"
  CHARACTER(LEN=filename_max) :: extpar_filename

  ! external parameter: dictionary which maps internal variable names
  ! onto GRIB2 shortnames or NetCDF var names.
  CHARACTER(LEN=filename_max) :: extpar_varnames_map_file

  ! information read from extpar file (no namelist parameter)
  !
  INTEGER ::  &           !< stores the landcover classification used for the external parameter data
    &  i_lctype(max_dom)  !< 1: Globcover2009, 2: GLC2000
                          !< defined in mo_ext_data_state:inquire_extpar_file

  INTEGER ::  &           !< number of landuse classes
    &  nclass_lu(max_dom)

  INTEGER ::  &           !< number of months in external data file
    &  nmonths_ext(max_dom)

  !!----------------------------------------------------------------------------

CONTAINS

  FUNCTION generate_filename(extpar_filename, model_base_dir, grid_filename, &
    &                        nroot, jlev, idom) &
    &  RESULT(result_str)
    CHARACTER(len=*), INTENT(IN)    :: extpar_filename, &
      &                                model_base_dir,  &
      &                                grid_filename
    INTEGER,          INTENT(IN)   :: nroot, jlev, idom
    CHARACTER(len=MAX_CHAR_LENGTH)  :: result_str
    TYPE (t_keyword_list), POINTER  :: keywords

    NULLIFY(keywords)
    CALL associate_keyword("<path>",     TRIM(model_base_dir), keywords)
    CALL associate_keyword("<gridfile>", TRIM(grid_filename),  keywords)
    CALL associate_keyword("<nroot>",  TRIM(int2string(nroot,"(i0)")),   keywords)
    CALL associate_keyword("<nroot0>", TRIM(int2string(nroot,"(i2.2)")), keywords)
    CALL associate_keyword("<jlev>",   TRIM(int2string(jlev, "(i2.2)")), keywords)
    CALL associate_keyword("<idom>",   TRIM(int2string(idom, "(i2.2)")), keywords)
    ! replace keywords in "extpar_filename", which is by default
    ! extpar_filename = "<path>extpar_<gridfile>"
    result_str = with_keywords(keywords, TRIM(extpar_filename))

  END FUNCTION generate_filename
!-----------------------------------------------------------------------
  FUNCTION generate_td_filename(extpar_td_filename, model_base_dir, grid_filename, month, year, clim) &
    &  RESULT(result_str)
    CHARACTER(len=*), INTENT(IN)    :: extpar_td_filename, &
      &                                model_base_dir,  &
      &                                grid_filename
    INTEGER, INTENT(IN)             :: month
    INTEGER, INTENT(IN), OPTIONAL   :: year
    LOGICAL, INTENT(IN), OPTIONAL   :: clim
    CHARACTER(len=4) :: syear
    CHARACTER(len=2) :: smonth
    CHARACTER(len=MAX_CHAR_LENGTH) :: result_str
    TYPE(t_keyword_list), POINTER :: keywords
    LOGICAL :: lclim
    CHARACTER(len=*), PARAMETER :: &
    &  routine = modname//':generate_td_filename:'

    lclim = .FALSE.
    IF (PRESENT(clim)) lclim = clim
    IF (PRESENT (year)) THEN
     WRITE(syear, '(i4.4)') year
    ELSEIF (lclim) THEN
      syear="CLIM"
    ELSE
      CALL finish(routine, 'Missing year for a non climatological run')
    END IF
    WRITE(smonth,'(i2.2)') month
    NULLIFY(keywords)
    CALL associate_keyword("<path>",     TRIM(model_base_dir), keywords)
    CALL associate_keyword("<gridfile>", TRIM(grid_filename),  keywords)
    CALL associate_keyword("<year>", syear,  keywords)
    CALL associate_keyword("<month>", smonth,  keywords)
    ! replace keywords in "extpar_filename", which is by default
    ! extpar_td_filename = "<path>extpar_<year>_<month>_<gridfile>"
    ! if clim ist present and clim=.TRUE., <year> ist subst. by "CLIM"
    result_str = TRIM(with_keywords(keywords, TRIM(extpar_td_filename)))

  END FUNCTION generate_td_filename
END MODULE mo_extpar_config
