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
  USE mo_exception,          ONLY: finish, message, message_text
  USE mo_cdi,                ONLY: cdi_undefid

  IMPLICIT NONE

  PRIVATE

  ! collection of external parameter file attributes
  ! for the atmosphere, which are read from the extpar file.
  TYPE :: t_ext_atm_attr
    INTEGER                    :: id              !< patch ID
    LOGICAL                    :: have_inquired   !< TRUE. if file content has been inquired
    CHARACTER(:), ALLOCATABLE  :: extpar_file     !< file name
    INTEGER                    :: cdi_extpar_id
    INTEGER                    :: cdi_filetype    !< CDI filetype (GRIB2/NetCDF)
    INTEGER                    :: nhori           !< number of horizon sectors
    INTEGER                    :: nclass_lu       !< number of landuse classes
    INTEGER                    :: nmonths_ext     !< number of months in external data file
    INTEGER                    :: i_lctype        !< landcover classification used for the
                                                  !  external parameter data
                                                  !  1: Globcover2009, 2: GLC2000
    LOGICAL                    :: is_frglac_in    !< TRUE, if the extpar file contains fr_glac
  CONTAINS
    ! initialize object with dummy values
    PROCEDURE  :: init         => ext_atm_attr__init
    ! print object values
    PROCEDURE  :: print_values => ext_atm_attr__print_values
  END TYPE t_ext_atm_attr

  ! collection of O3 file attributes, which are read from the o3 file.
  TYPE :: t_ext_o3_attr
    INTEGER                    :: id              !< patch ID
    LOGICAL                    :: have_inquired   !< TRUE. if file content has been inquired
    CHARACTER(:), ALLOCATABLE  :: levelname
    CHARACTER(:), ALLOCATABLE  :: cellname
    CHARACTER(:), ALLOCATABLE  :: o3name
    CHARACTER(:), ALLOCATABLE  :: o3unit
    INTEGER                    :: nlev_o3         !< number of vertical levels
    INTEGER                    :: nmonths         !< number of months in O3 file
  CONTAINS
    ! initialize object with dummy values
    PROCEDURE  :: init         => ext_o3_attr__init
    ! print object values
    PROCEDURE  :: print_values => ext_o3_attr__print_values
  END TYPE t_ext_o3_attr

  ! Number of landcover classes provided by external parameter data
  ! Needs to be changed into a variable if landcover classifications
  ! with a different number of classes become available
  INTEGER, PARAMETER :: num_lcc = 23, n_param_lcc = 8

  ! types
  PUBLIC :: t_ext_atm_attr
  PUBLIC :: t_ext_o3_attr

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
  PUBLIC :: ext_atm_attr
  PUBLIC :: ext_o3_attr
  PUBLIC :: num_lcc
  PUBLIC :: n_param_lcc

  ! subroutines/functions
  PUBLIC :: generate_filename
  PUBLIC :: generate_td_filename

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_extpar_config'


  !>
  !!----------------------------------------------------------------------------
  !! Derived type containing control variables specific to the nonhydrostatic
  !! atm model
  !!------------------------------------------------------------------------

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

  ! ExtPar input filename, may contain keywords, by default
  ! extpar_filename = "<path>extpar_<gridfile>"
  CHARACTER(LEN=filename_max) :: extpar_filename

  ! external parameter: dictionary which maps internal variable names
  ! onto GRIB2 shortnames or NetCDF var names.
  CHARACTER(LEN=filename_max) :: extpar_varnames_map_file

  TYPE(t_ext_atm_attr), TARGET:: ext_atm_attr(max_dom)  !< file attributes read from extpar file
  TYPE(t_ext_o3_attr),  TARGET:: ext_o3_attr(max_dom)   !< file attributes read from o3 file
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


  !>
  !! Initializes object of type t_ext_atm_attr
  !!
  !! Initializes object of type t_ext_atm_attr
  !! with dummy values.
  !!
  SUBROUTINE ext_atm_attr__init (this, patch_id)
    CLASS(t_ext_atm_attr) :: this
    INTEGER, INTENT(IN) :: patch_id

    this%id              = patch_id
    this%have_inquired   = .FALSE.
    this%extpar_file     = ""
    this%cdi_extpar_id   = cdi_undefid
    this%cdi_filetype    = cdi_undefid
    this%nhori           = -1
    this%nclass_lu       = -1
    this%nmonths_ext     = -1
    this%i_lctype        = -1
    this%is_frglac_in    = .FALSE.
  END SUBROUTINE ext_atm_attr__init


  !>
  !! Print object of type t_ext_atm_attr
  !!
  !! Print values for object of type t_ext_atm_attr
  !!
  SUBROUTINE ext_atm_attr__print_values (this)
    CLASS(t_ext_atm_attr) :: this

    WRITE(message_text, '(a,i2,a)'   ) "ext_atm_attr for DOM ", this%id, ":"
    CALL message ('  ',message_text)
    !
    WRITE(message_text, '(a,l1)' ) "have_inquired: " , this%have_inquired
    CALL message ('  ',message_text)
    WRITE(message_text, '(a,a)' ) "extpar_file: "    , this%extpar_file
    CALL message ('  ',message_text)
    WRITE(message_text, '(a,i2)') "cdi_extpar_id: "  , this%cdi_extpar_id
    CALL message ('  ',message_text)
    WRITE(message_text, '(a,i2)') "cdi_filetype: "   , this%cdi_filetype
    CALL message ('  ',message_text)
    WRITE(message_text, '(a,i4)') "nhori: "          , this%nhori
    CALL message ('  ',message_text)
    WRITE(message_text, '(a,i2)') "nclass_lu: "      , this%nclass_lu
    CALL message ('  ',message_text)
    WRITE(message_text, '(a,i2)') "nmonths_ext: "    , this%nmonths_ext
    CALL message ('  ',message_text)
    WRITE(message_text, '(a,i2)') "i_lctype: "       , this%i_lctype
    CALL message ('  ',message_text)
    WRITE(message_text, '(a,l1)' ) "is_frglac_in: "  , this%is_frglac_in
    CALL message ('  ',message_text)
    CALL message ('  ','')

  END SUBROUTINE ext_atm_attr__print_values


  !>
  !! Initializes object of type t_ext_o3_attr
  !!
  !! Initializes object of type t_ext_o3_attr
  !! with dummy values.
  !!
  SUBROUTINE ext_o3_attr__init (this, patch_id)
    CLASS(t_ext_o3_attr) :: this
    INTEGER, INTENT(IN) :: patch_id

    this%id              = patch_id
    this%have_inquired   = .FALSE.
    this%levelname       = ""
    this%cellname        = ""
    this%o3name          = ""
    this%o3unit          = ""
    this%nlev_o3         = -1
    this%nmonths         = -1
  END SUBROUTINE ext_o3_attr__init

  !>
  !! Print object of type t_ext_o3_attr
  !!
  !! Print values for object of type t_ext_o3_attr
  !!
  SUBROUTINE ext_o3_attr__print_values (this)
    CLASS(t_ext_o3_attr) :: this

    WRITE(message_text, '(a,i2,a)'   ) "ext_o3_attr for DOM ", this%id, ":"
    CALL message ('  ',message_text)
    !
    WRITE(message_text, '(a,l1)' ) "have_inquired: ", this%have_inquired
    CALL message ('  ',message_text)
    WRITE(message_text, '(a,a)' ) "levelname: "     , this%levelname
    CALL message ('  ',message_text)
    WRITE(message_text, '(a,a)') "cellname: "       , this%cellname
    CALL message ('  ',message_text)
    WRITE(message_text, '(a,a)') "o3name: "         , this%o3name
    CALL message ('  ',message_text)
    WRITE(message_text, '(a,a)') "o3unit: "         , this%o3unit
    CALL message ('  ',message_text)
    WRITE(message_text, '(a,i2)') "nlev_o3: "       , this%nlev_o3
    CALL message ('  ',message_text)
    WRITE(message_text, '(a,i2)') "nmonths: "       , this%nmonths
    CALL message ('  ',message_text)
    CALL message ('  ','')

  END SUBROUTINE ext_o3_attr__print_values

END MODULE mo_extpar_config
