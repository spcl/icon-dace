! NWP surface tile utility functions
!
! Surface tile utility functions. Many of them are dealing with the
! request of tile metainformation and its conversion between the
! model-internal tile structure and the output tile structure according
! to GRIB2 Product Definition Template (PDT) 4.55.
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

MODULE mo_nwp_sfc_tiles

  USE mo_impl_constants,     ONLY: MAX_CHAR_LENGTH
  USE mo_exception,          ONLY: finish, message, message_text
  USE mo_util_string,        ONLY: int2string
  USE mo_util_table,         ONLY: t_table, initialize_table, add_table_column, &
    &                              set_table_entry, print_table, finalize_table
  USE mo_mpi,                ONLY: my_process_is_stdio
  USE mo_var_metadata_types, ONLY: CLASS_TILE, CLASS_TILE_LAND


  IMPLICIT NONE

  PRIVATE

  ! types
  PUBLIC :: t_tile_list
  PUBLIC :: t_tile_att
  PUBLIC :: t_tileinfo_icon
  PUBLIC :: t_tileinfo_grb2

  ! subroutines
  PUBLIC :: setup_tile_list

  ! variables
  PUBLIC :: trivial_tile_att


  !> module name string
  CHARACTER(LEN = *), PARAMETER :: modname = 'mo_nwp_sfc_tiles'


  !>
  !! tile attributes 
  !!
  TYPE :: t_tile_att
    INTEGER :: tile_id_icon       !< ICON internal tile ID
    INTEGER :: tile_att           !< attribute according to GRIB2 code table 4.241
    !
    CLASS(t_tile), POINTER :: tile_ptr !< concatenates attributes with its 
                                       !< corresponding tile

  CONTAINS
    PROCEDURE  :: getTileinfo_icon => tile_att_getTileinfo_icon
    PROCEDURE  :: getTileinfo_grb2 => tile_att_getTileinfo_grb2

    ! convert tile/attribute pair into ICON-internal varname suffix
    PROCEDURE  :: getTileSuffix => tile_att_getTileSuffix
  END TYPE t_tile_att


  !>
  !! individual tile meta-information 
  !!
  TYPE :: t_tile
    CHARACTER(len=MAX_CHAR_LENGTH) :: name      !< tile name
    INTEGER                        :: tile_id
    TYPE(t_tile_att), ALLOCATABLE  :: atts(:)   !< array of tile attributes  

  CONTAINS
    !
    ! initialize new tile
    PROCEDURE  :: init                      => tile_init
    !
    ! append tile attribute
    PROCEDURE  :: append_att                => tile_append_att
    !
    ! get total number of tile attributes 
    PROCEDURE  :: getNumberOfTileAttributes => tile_getNumberOfTileAttributes
    !
    ! finalization
    PROCEDURE  :: destruct                  => tile_destruct
  END TYPE t_tile


  !>
  !! list of tiles
  !!
  TYPE :: t_tile_list
    INTEGER                   :: ntiles_lnd   !< number of land tiles
    TYPE(t_tile), ALLOCATABLE :: tile(:)      !< tiles

  CONTAINS
    !
    ! get number of tiles
    PROCEDURE  :: getNumberOfTiles => tile_list_getNumberOfTiles

    ! given grib2/icon tile info, return corresponding tile attribute object
    PROCEDURE  :: getTileAtt       => tile_list_getTileAtt

    ! visualize basic surface tile setup
    PROCEDURE  :: printSetup       => tile_list_printSetup

    ! wrapper routines for tile-specific requests
    PROCEDURE  :: getTileinfo_icon          => tile_list_getTileinfo_icon
    PROCEDURE  :: getTileinfo_grb2          => tile_list_getTileinfo_grb2
    PROCEDURE  :: getNumberOfTileAttributes => tile_list_getNumberOfTileAttributes

    ! finalization
    PROCEDURE  :: destruct         => tile_list_destruct
  END TYPE t_tile_list


  !> tile information, abstract base class
  TYPE, ABSTRACT :: t_tileinfo_elt
  END TYPE t_tileinfo_elt

  !> tile information, ICON-internal indexing
  !!
  TYPE, EXTENDS(t_tileinfo_elt) :: t_tileinfo_icon
    INTEGER :: idx     !< tile index (variable specific)
  END TYPE t_tileinfo_icon

  !> tile information, according to GRIB2 Product Definition Template (PDT) 4.55
  !!
  TYPE, EXTENDS(t_tileinfo_elt) :: t_tileinfo_grb2
    INTEGER :: idx     !< tile index (variable specific)
    INTEGER :: att     !< tile attribute (variable specific)
  END TYPE t_tileinfo_grb2


  ! trivial tile information, denoting a "no-tile" field:
  TYPE(t_tile),          TARGET    :: trivial_tile
  TYPE(t_tile_att),      TARGET    :: trivial_tile_att

  TYPE(t_tileinfo_grb2), PARAMETER :: trivial_tileinfo_grb2  = t_tileinfo_grb2(idx = 0, att = 0)
  INTEGER,               PARAMETER :: trivial_tileId         = 0

  !
  ! maximum number of tile attributes
  INTEGER, PARAMETER :: max_atts = 4

CONTAINS


  !! Setup list of active surface tiles
  SUBROUTINE setup_tile_list (tile_list, ntiles_lnd, lsnowtile, isub_water, isub_lake, isub_seaice)

    TYPE(t_tile_list), TARGET, INTENT(INOUT) :: tile_list !< list of tiles

    INTEGER,    INTENT(IN) :: ntiles_lnd       !< number of land tiles
    LOGICAL,    INTENT(IN) :: lsnowtile        !< snowtiles yes/no

    INTEGER,    INTENT(IN) :: isub_water       !< internal array positions for water
    INTEGER,    INTENT(IN) :: isub_lake        !< lake
    INTEGER,    INTENT(IN) :: isub_seaice      !< seaice


    ! local
    CHARACTER(LEN=*), PARAMETER :: routine = modname//":setup_tile_list"
    INTEGER :: istat            ! error flag
    INTEGER :: ntiles           ! total number of spatio-temporal changing tiles
    INTEGER :: i                ! loop index

    ! coverage attributes according to GRIB2 code table 4.241
    INTEGER, PARAMETER :: UNDEF = 0  ! undefined
    INTEGER, PARAMETER :: UNMOD = 1  ! unmodified
    INTEGER, PARAMETER :: SNOW  = 2  ! snow-covered
    INTEGER, PARAMETER :: SEAICE= 4  ! sea ice covered
    !-----------------------------------------------------------------------

    ! first, initialize "trivial" tile:
    CALL trivial_tile%init(tile_id=trivial_tileinfo_grb2%idx)
    CALL trivial_tile%append_att(tile_id_icon=trivial_tileId, tile_att=trivial_tileinfo_grb2%att)
    !
    ! initialize trivial_tile_att
    trivial_tile_att = trivial_tile%atts(1)
    !
    ! Paranoia
    IF (.NOT.ASSOCIATED(trivial_tile_att%tile_ptr)) THEN
      CALL finish(routine, 'Initialization of trivial_tile_att failed. trivial_tile_att%tile_ptr not associated')
    ENDIF

    !
    ! ASCII art for ICON-internal tile structure, depicting the configuration for:
    ! - ntiles_lnd = 3
    ! - lsnowtile  = .TRUE.
    !
    !------------------------------------------------------------------!
    !   1  |   2  |   3  ||   4  |   5  |   6  ||   7  |   8  |    9   !
    ! land | land | land || snow | snow | snow || open | lake | seaice !
    !      |      |      ||      |      |      || sea  |      |        !
    !------------------------------------------------------------------!
    !
    !<---  ntiles_lnd --->                      <---- ntiles_water ---->
    !
    !<-------------- ntiles_total ------------>
    !
    !<---------------- ntiles_total + ntiles_water ------------------->
    !

    ! ASCII art for corresponding GRIB2 tile structure:
    !
    !-----------------------------------------------------------------!
    !    land    |    land    |    land    |     oce    |    lake     !
    !      1     |      2     |      3     |      4     |     5       !  <- GRIB2 Tile ID
    !     / \    |     / \    |     / \    |     / \    |     |       !
    !****/***\***|****/***\***|****/***\***|****/***\***|*****|*******!  *******************
    ! umod | snw | umod | snw | umod | snw | umod | ice |   undef     !  <- Attribute
    !  1   |  4  |   2  |  5  |   3  |  6  |   7  |  9  |     8       !  <- Internal Tile ID
    !-----------------------------------------------------------------!     


    IF (ntiles_lnd>1) THEN
      ! The total number of spatio-temporal changing tiles (see also GRIB2 PDT 4.55) 
      ! of an individual grid cell is given as: 
      ! the total number of land tiles + one ocean tile + one lake tile  
      ntiles = ntiles_lnd + 2
    ELSE
      ntiles = ntiles_lnd
    ENDIF
    !
    ! allocate list of tiles
    ALLOCATE(tile_list%tile(ntiles), STAT=istat)
    IF (istat /= 0) THEN
      CALL finish(routine, &
        &      'allocation of tile_list failed')
    ENDIF

    IF (ntiles>1) THEN

      ! store number of land tiles
      tile_list%ntiles_lnd = ntiles_lnd

      ! initialize tiles
      DO i=1,ntiles
        CALL tile_list%tile(i)%init(tile_id=i)
      ENDDO

      ! append tile attributes
      !
      ! we adhere to the ASCII art above and start with the land tiles
      DO i=1,ntiles_lnd
        CALL tile_list%tile(i)%append_att(tile_id_icon=i, tile_att=MERGE(UNMOD, UNDEF, lsnowtile))
        IF (lsnowtile) THEN
          ! If snowtiles exist, a second set of attributes is added (snow covered)
          CALL tile_list%tile(i)%append_att(tile_id_icon=i+ntiles_lnd, tile_att=SNOW)
        ENDIF
        ! add name for beautification
        tile_list%tile(i)%name=TRIM("Land")
      ENDDO
      !
      ! tile with tile_id=ntiles_lnd+1 is defined as ocean tile
      CALL tile_list%tile(ntiles_lnd+1)%append_att(tile_id_icon=isub_water, tile_att=UNMOD)
      CALL tile_list%tile(ntiles_lnd+1)%append_att(tile_id_icon=isub_seaice, tile_att=SEAICE)
      ! add name for beautification
      tile_list%tile(ntiles_lnd+1)%name=TRIM("Ocean")
      !
      ! tile with tile_id=ntiles_lnd+2 is defined as lake tile
      CALL tile_list%tile(ntiles_lnd+2)%append_att(tile_id_icon=isub_lake, tile_att=UNDEF)
      ! add name for beautification
      tile_list%tile(ntiles_lnd+2)%name=TRIM("Lake")
    ELSE

      ! 'trivial' tile

      ! store dummy value for number of land tiles
      tile_list%ntiles_lnd = 0

      ! initialize 'trivial' tile
      CALL tile_list%tile(1)%init(tile_id=trivial_tileinfo_grb2%idx)
      CALL tile_list%tile(1)%append_att(tile_id_icon=1, tile_att=trivial_tileinfo_grb2%att)
    ENDIF

    ! debug printput
    CALL tile_list%printSetup()

  END SUBROUTINE setup_tile_list


  !> return the corresponding ICON tile info
  !!
  FUNCTION tile_att_getTileinfo_icon(tile_att) RESULT(tileinfo_icon)
    TYPE(t_tileinfo_icon) :: tileinfo_icon
    CLASS(t_tile_att), TARGET, INTENT(IN) :: tile_att  !< passed-object dummy argument

    tileinfo_icon%idx = tile_att%tile_id_icon
  END FUNCTION tile_att_getTileinfo_icon


  !> return the corresponding GRIB2 tile info
  !!
  FUNCTION tile_att_getTileinfo_grb2(tile_att) RESULT(tileinfo_grb2)
    TYPE(t_tileinfo_grb2) :: tileinfo_grb2
    CLASS(t_tile_att), TARGET, INTENT(IN) :: tile_att  !< passed-object dummy argument

    tileinfo_grb2%idx  = tile_att%tile_ptr%tile_id
    tileinfo_grb2%att  = tile_att%tile_att
  END FUNCTION tile_att_getTileinfo_grb2


  !> Given a GRIB2/ICON tileinfo object, return the corresponding ICON tile info (wrapper routine)
  !!
  FUNCTION tile_list_getTileinfo_icon(tile_list, tileinfo) RESULT(tileinfo_icon)
    TYPE(t_tileinfo_icon) :: tileinfo_icon
    CLASS(t_tile_list), TARGET, INTENT(IN) :: tile_list  !< passed-object dummy argument
    CLASS(t_tileinfo_elt), INTENT(IN) :: tileinfo
    ! local
    TYPE(t_tile_att), POINTER   :: tile_att

    tile_att => tile_list%getTileAtt(tileinfo)
    tileinfo_icon = tile_att%getTileinfo_icon()
  END FUNCTION tile_list_getTileinfo_icon


  !> Given a GRIB2/ICON tileinfo object, return the corresponding GRB2 tile info (wrapper routine)
  !!
  FUNCTION tile_list_getTileinfo_grb2(tile_list, tileinfo) RESULT(tileinfo_grb2)
    TYPE(t_tileinfo_grb2) :: tileinfo_grb2
    CLASS(t_tile_list), TARGET, INTENT(IN) :: tile_list  !< passed-object dummy argument
    CLASS(t_tileinfo_elt), INTENT(IN) :: tileinfo
    ! local
    TYPE(t_tile_att), POINTER   :: tile_att

    tile_att => tile_list%getTileAtt(tileinfo)
    tileinfo_grb2 = tile_att%getTileinfo_grb2()
  END FUNCTION tile_list_getTileinfo_grb2




  !! For a given tile attribute, return ICON-internal varname suffix '_t_X'
  !! For a trivial tile attribute it returns an empty string.
  !!
  FUNCTION tile_att_getTileSuffix(tile_att) RESULT (tileSuffix)

    CLASS(t_tile_att),     INTENT(IN) :: tile_att  !< passed-object dummy argument

    TYPE(t_tileinfo_icon)      :: tileinfo_icon
    INTEGER                    :: icon_id
    CHARACTER(LEN=2)           :: icon_id_str
    CHARACTER(LEN=5)           :: tileSuffix
    !-----------------------------------------------------------------------

    ! get the tile ID
    tileinfo_icon = tile_att%getTileinfo_icon()
    icon_id       = tileinfo_icon%idx

    IF (icon_id == trivial_tileId) THEN
      ! special case: trivial tileID => empty suffix
      tileSuffix =''
    ELSE
      ! convert tile ID to the corresponding suffix string
      WRITE(icon_id_str,'(i2)') icon_id
      tileSuffix = '_t_'//TRIM(ADJUSTL(icon_id_str))
    ENDIF
  END FUNCTION tile_att_getTileSuffix

  
  !! Initializes tile of type t_tile
  !!
  SUBROUTINE tile_init (tile, tile_id, opt_name)

    CLASS(t_tile), INTENT(INOUT) :: tile     !< passed-object dummy argument

    INTEGER,       INTENT(IN)    :: tile_id
    CHARACTER(len=MAX_CHAR_LENGTH), OPTIONAL, INTENT(IN) :: opt_name

    ! local
    CHARACTER(LEN=*), PARAMETER :: routine = modname//":tile_init"
    INTEGER :: istat             ! error flag
    INTEGER :: iatt              ! loop index
  !-----------------------------------------------------------------------

    ! allocate array of tile attributes
    ALLOCATE(tile%atts(max_atts), STAT=istat)
    IF (istat /= 0) THEN
      CALL finish(routine, 'allocation of atts failed')
    ENDIF

    ! init name
    IF (PRESENT(opt_name)) THEN
      tile%name = TRIM(opt_name)
    ELSE
      tile%name = "undef"
    ENDIF

    ! store tile_id
    tile%tile_id = tile_id

    ! nullify pointer in order to have a defined association status
    DO iatt=1,SIZE(tile%atts)
      tile%atts(iatt)%tile_ptr => NULL()
    ENDDO

  END SUBROUTINE tile_init


  !! Appends set of attributes to tile. In addition 
  !! the new attribute set is concatenated with the 'parent' 
  !! tile by pointing to it.
  !!
  SUBROUTINE tile_append_att (tile, tile_id_icon, tile_att)

    CLASS(t_tile), TARGET, INTENT(INOUT) :: tile    !< passed-object dummy argument

    INTEGER,       INTENT(IN)    :: tile_id_icon    !< ICON internal tile ID
    INTEGER,       INTENT(IN)    :: tile_att        !< tile attribute 
                                                    ! see GRIB2 code table 4.241

    ! local
    INTEGER :: istat             ! error flag
    INTEGER :: iatt              ! loop index
    CHARACTER(LEN=*), PARAMETER :: routine = modname//":tile_append_att"
  !-----------------------------------------------------------------------

    istat = -1

    ! find next free element in attribute array and store attributes
    !
    DO iatt=1,SIZE(tile%atts)

      IF (ASSOCIATED(tile%atts(iatt)%tile_ptr)) THEN
        ! element atts(iatt) already in use
        CYCLE
      ELSE
        ! store metadata
        tile%atts(iatt)%tile_id_icon = tile_id_icon
        tile%atts(iatt)%tile_att     = tile_att

        ! point to parent tile
        tile%atts(iatt)%tile_ptr => tile

        istat = 0

        EXIT
      ENDIF
    ENDDO

    IF (istat /= 0) THEN
      CALL finish(routine, 'appending tile attributes failed')
    ENDIF
  END SUBROUTINE tile_append_att



  !! Get total number of attributes for a given tile
  !!
  INTEGER FUNCTION tile_getNumberOfTileAttributes(tile) RESULT(natt)

    CLASS(t_tile), INTENT(IN) :: tile     !< passed-object dummy argument

    ! local
    INTEGER :: iatt   ! loop index
  !-----------------------------------------------------------------------

    natt = 0
    DO iatt=1,SIZE(tile%atts)
      IF (ASSOCIATED(tile%atts(iatt)%tile_ptr)) THEN
        natt = natt + 1
      ENDIF
    ENDDO

  END FUNCTION tile_getNumberOfTileAttributes



  !> Given the GRIB2/ICON tileinfo, return the number of tile attributes (wrapper routine)
  !!
  FUNCTION tile_list_getNumberOfTileAttributes(tile_list, tileinfo) RESULT(natt)
    INTEGER :: natt
    CLASS(t_tile_list), TARGET, INTENT(IN) :: tile_list  !< passed-object dummy argument
    CLASS(t_tileinfo_elt), INTENT(IN) :: tileinfo
    ! local
    TYPE(t_tile_att), POINTER   :: this_att
    TYPE(t_tile),     POINTER   :: this_tile

    this_att  => tile_list%getTileAtt(tileinfo)
    this_tile => this_att%tile_ptr
    natt      =  this_tile%getNumberOfTileAttributes()
  END FUNCTION tile_list_getNumberOfTileAttributes


  !! Destructor for variable of type t_tile
  !!
  SUBROUTINE tile_destruct (tile)
    CLASS(t_tile), INTENT(INOUT) :: tile     !< passed-object dummy argument

    ! local
    INTEGER :: istat            ! error flag
    CHARACTER(LEN=*), PARAMETER :: routine = modname//":tile_destruct"
  !-----------------------------------------------------------------------

    DEALLOCATE(tile%atts, STAT=istat)
    IF (istat /= 0) THEN
      CALL finish(routine, 'deallocation of tile failed')
    ENDIF

  END SUBROUTINE tile_destruct


  !>
  !! Get total number of tiles
  !!
  !! Get variable-specific number of tiles. All tiled fields 
  !! belong to one of the following classes: 
  !! CLASS_TILE: variable contains land and water tiles
  !! CLASS_TILE_LAND: variable contains only land tiles
  !! This information is part of each variables metadata (see 'info').
  !!
  !! Note that the number of tiles returned by this routine adheres to the 
  !! output (GRIB2) tile structure and not the model-internal structure. 
  !! I.e. snowtiles and the sea-ice tile are not considered as separate tiles.
  !!
  INTEGER FUNCTION tile_list_getNumberOfTiles(tile_list,varClass) RESULT(numberOfTiles)

    CLASS(t_tile_list), INTENT(IN) :: tile_list  !< passed-object dummy argument

    INTEGER,            INTENT(IN) :: varClass   !< variable class
                                                 !< CLASS_TILE
                                                 !<   variable contains land and water tiles 
                                                 !< CLASS_TILE_LAND
                                                 !<   variable contains only land tiles

    ! local
    CHARACTER(LEN=*), PARAMETER :: routine = modname//":tile_list_getNumberOfTiles"
  !-----------------------------------------------------------------------

    SELECT CASE(varClass)
    CASE(CLASS_TILE)
      numberOfTiles = SIZE(tile_list%tile)
    CASE(CLASS_TILE_LAND)
      numberOfTiles = tile_list%ntiles_lnd
    CASE DEFAULT
      CALL finish( routine, 'class not supported' )
    END SELECT

  END FUNCTION tile_list_getNumberOfTiles


  !>
  !! Given a tile info object, return the corresponding tile attribute object.
  !!
  FUNCTION tile_list_getTileAtt(tile_list, tileinfo) RESULT(this_att)
    TYPE(t_tile_att), POINTER   :: this_att  ! pointer to attribute

    CLASS(t_tile_list), TARGET, INTENT(IN) :: tile_list  !< passed-object dummy argument
    CLASS(t_tileinfo_elt),      INTENT(IN) :: tileinfo

    ! local
    CHARACTER(LEN=*), PARAMETER :: routine = modname//":tile_list_getTileAtt"
    INTEGER                     :: i, icon_id, iatt
    TYPE(t_tile),     POINTER   :: this_tile

    LOGICAL                     :: lfound    ! matching tile info found (TRUE/FALSE)

    !-----------------------------------------------------------------------

    SELECT TYPE(tileinfo)
    TYPE IS (t_tileinfo_icon)
      icon_id = tileinfo%idx
      IF (icon_id == trivial_tileId) THEN
        this_att => trivial_tile_att
      ELSE IF (icon_id < 0) THEN
        CALL finish(routine, "illegal ICON tileID "//int2string(icon_id))
      ELSE
        lfound=.FALSE.
        !
        loop_tile1: DO i=1, SIZE(tile_list%tile)
          this_tile => tile_list%tile(i)
          !
          loop_att1: DO iatt=1, SIZE(this_tile%atts)
            this_att => this_tile%atts(iatt)
            IF (.NOT.ASSOCIATED(this_att%tile_ptr)) CYCLE
            !
            IF (this_att%tile_id_icon == icon_id) THEN
              lfound=.TRUE.
              EXIT loop_tile1
            ENDIF
          ENDDO loop_att1
        ENDDO loop_tile1
        !
        IF (.NOT.lfound) THEN
          CALL finish(routine, "no list entry found for ICON tileID "//int2string(icon_id))
        ENDIF
      ENDIF
      !
    TYPE IS (t_tileinfo_grb2)
      lfound=.FALSE.
      !
      loop_tile2: DO i=1, SIZE(tile_list%tile)
        this_tile => tile_list%tile(i)
        !
        loop_att2: DO iatt=1, SIZE(this_tile%atts)
          this_att => this_tile%atts(iatt)
          IF (.NOT.ASSOCIATED(this_att%tile_ptr)) CYCLE
          !
          IF ((this_att%tile_ptr%tile_id == tileinfo%idx) .AND. &
            & (this_att%tile_att         == tileinfo%att)) THEN
            lfound=.TRUE.
            EXIT loop_tile2
          ENDIF
        ENDDO loop_att2
      ENDDO loop_tile2
      !
      IF (.NOT.lfound) THEN
        CALL finish(routine, "no list entry found for GRIB2 tile "//&
          &int2string(tileinfo%idx)//", "//int2string(tileinfo%att))
      ENDIF
    END SELECT
  END FUNCTION tile_list_getTileAtt


  !! Visualize basic surface tile setup. 
  !! Graphical visualization of mapping between internal tile ids 
  !! and grib2 tile information.
  !!
  SUBROUTINE tile_list_printSetup (tile_list)
    !
    CLASS(t_tile_list), TARGET, INTENT(IN) :: tile_list   !< passed-object dummy argument

    ! local variables
    INTEGER         :: it              ! loop counter
    INTEGER         :: iatt            ! loop counter
    TYPE(t_table)   :: table
    INTEGER         :: irow            ! row to fill
    !
    CHARACTER(LEN = *), PARAMETER :: tileNameCol   = "tileName", &
      &                              tileIdCol     = "GRIB2 tileId"  , &
      &                              tileAttCol    = "tileAtts", &
      &                              tileIdIconCol = "ICON tile Ids"


    TYPE(t_tile),     POINTER :: this_tile
    CHARACTER(len=MAX_CHAR_LENGTH) :: tileName

    CHARACTER(len=MAX_CHAR_LENGTH) :: att_str  ! attribute string
    CHARACTER(len=MAX_CHAR_LENGTH) :: id_str   ! id string

    CHARACTER(LEN=*), PARAMETER :: routine = modname//":tile_list_printSetup"
    !--------------------------------------------------------------------------

    ! will only be executed by stdio process
    IF(.NOT. my_process_is_stdio()) RETURN

    ! poor man's table header 
    WRITE(message_text,'(a,a,a,i2)') 'Surface tile setup'
    CALL message('', TRIM(message_text))


    ! table-based output
    CALL initialize_table(table)
    ! the latter is no longer mandatory
    CALL add_table_column(table, tileNameCol)
    CALL add_table_column(table, tileIdCol)
    CALL add_table_column(table, tileAttCol)
    CALL add_table_column(table, tileIdIconCol)



    ! initialize counter
    irow = 0

    DO it=1,SIZE(tile_list%tile)
      ! print tile setup
      this_tile => tile_list%tile(it)

      irow = irow + 1
      !
      ! print name
      tileName = this_tile%name 
      CALL set_table_entry(table,irow,tileNameCol, ADJUSTL(TRIM(tileName)))
      !
      ! print tile id
      CALL set_table_entry(table,irow,tileIdCol, TRIM(int2string(this_tile%tile_id)))
      !
      ! print attributes
      att_str=""
      DO iatt=1,SIZE(this_tile%atts)
       IF (.NOT.ASSOCIATED(this_tile%atts(iatt)%tile_ptr)) CYCLE
       att_str = TRIM(att_str)//TRIM(int2string(this_tile%atts(iatt)%tile_att))//","
      ENDDO
      CALL set_table_entry(table,irow,tileAttCol, ADJUSTL(TRIM(att_str)))
      !
      ! print internal ICON id
      id_str=""
      DO iatt=1,SIZE(this_tile%atts)
       IF (.NOT.ASSOCIATED(this_tile%atts(iatt)%tile_ptr)) CYCLE
       id_str = TRIM(id_str)//TRIM(int2string(this_tile%atts(iatt)%tile_id_icon))//","
      ENDDO
      CALL set_table_entry(table,irow,tileIdIconCol, ADJUSTL(TRIM(id_str)))

    ENDDO

    CALL print_table(table, opt_delimiter=' | ')
    CALL finalize_table(table)

    WRITE (0,*) " " ! newline
  END SUBROUTINE tile_list_printSetup


  !! Destructor for variable of type t_tile_list
  !! Note that we have voted against a type-bound procedure 
  !! and binding name here. Likely, this makes it easier to 
  !! identify the position in the code where this routine is 
  !! called.
  !!
  SUBROUTINE tile_list_destruct (tile_list)

    CLASS(t_tile_list), INTENT(INOUT) :: tile_list     !< passed-object dummy argument

    ! local
    INTEGER :: istat            ! error flag
    INTEGER :: i                ! loop index
    CHARACTER(LEN=*), PARAMETER :: routine = modname//":tile_list_destruct"
    !-----------------------------------------------------------------------

    DO i=1,SIZE(tile_list%tile)
      CALL tile_list%tile(i)%destruct()
    ENDDO

    DEALLOCATE(tile_list%tile, STAT=istat)
    IF (istat /= 0) THEN
      CALL finish(routine, 'deallocation of tile_list failed')
    ENDIF

  END SUBROUTINE tile_list_destruct

END MODULE mo_nwp_sfc_tiles

