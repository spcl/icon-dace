!> Contains methods for initialization of JSBACH usecases.
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2024, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
MODULE mo_jsb_model_usecases
#ifndef __NO_JSBACH__

  USE mo_exception,           ONLY: message, finish

  USE mo_jsb_model_class,     ONLY: t_jsb_model
  USE mo_jsb_process_class,   ONLY: A2L_, L2A_, SEB_, TURB_, SSE_, HYDRO_, RAD_, HD_, ASSIMI_, PHENO_, CARBON_, DISTURB_, &
                                    FUEL_, ALCC_, NLCC_, VEG_, PPLCC_, TLCC_, TCQ_
  USE mo_jsb_process_class,   ONLY: ON_LEAFS_, ON_TILE_, AGGREGATE_ !, ON_SUBTREE_
  USE mo_jsb_tile,            ONLY: t_jsb_tile
  USE mo_jsb_lct_class,       ONLY: t_jsb_lct, &
    &                               LAND_TYPE, VEG_TYPE, GLACIER_TYPE, LAKE_TYPE !, BARE_TYPE

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: init_usecase

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_model_usecases'

CONTAINS

  SUBROUTINE init_usecase(model)

    TYPE(t_jsb_model), POINTER, INTENT(inout) :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':init_usecase'

    IF (model%config%use_lakes) THEN
      CALL message(TRIM(routine), 'Initializing JSBACH usecase "'//TRIM(model%config%usecase)//'" with lakes')
    ELSE
      CALL message(TRIM(routine), 'Initializing JSBACH usecase "'//TRIM(model%config%usecase)//'"')
    END IF


    SELECT CASE (TRIM(model%config%usecase))

    CASE ('jsbach_lite')
      IF (model%config%use_tmx) THEN
        CALL init_usecase_lite_new(model)
      ELSE
        CALL init_usecase_lite(model)
      END IF

    CASE ('jsbach_pfts')
      CALL init_usecase_pfts(model)

    CASE DEFAULT
      CALL finish(TRIM(routine), 'JSBACH usecase not defined.')

    END SELECT

  END SUBROUTINE init_usecase

  SUBROUTINE init_usecase_lite(model)

    TYPE(t_jsb_model), POINTER, INTENT(inout) :: model

    CLASS(t_jsb_tile), POINTER :: box_tile, lake_tile, land_tile, glacier_tile, veg_tile

    CHARACTER(len=*), PARAMETER :: routine = modname//':init_usecase_lite'

    ! Create root of the tile structure representing the whole grid box
    box_tile => t_jsb_tile('box', 'Top tile', &
      & processes      =(/A2L_,     L2A_,     RAD_,      HYDRO_,    TURB_    /), &
      & process_actions=(/ON_TILE_, ON_TILE_, ON_LEAFS_, ON_LEAFS_, ON_LEAFS_/), &
      & model_id=model%id, fract_varname='notsea')
    model%top => box_tile

#ifndef __NO_JSBACH_HD__
    IF (model%processes(HD_)%p%config%active) THEN
      CALL box_tile%Set_process_action(HD_, ON_TILE_)
    END IF
#endif

      ! Create sub-tiles
      IF (model%config%use_lakes) THEN
        lake_tile => t_jsb_tile('lake', 'Lake tile', &
          & processes=(/SEB_/), process_actions=(/ON_TILE_/), &
          & lct=t_jsb_lct(LAKE_TYPE), parent=box_tile, fract_varname='fract_lake')
        land_tile => t_jsb_tile('land', 'Land tile', &
          & processes=(/SEB_/), process_actions=(/ON_TILE_/), &
          & lct=t_jsb_lct(LAND_TYPE), parent=box_tile)
      ELSE
        land_tile => t_jsb_tile('land', 'Land tile', &
          & processes=(/SEB_/), process_actions=(/ON_TILE_/), &
          & lct=t_jsb_lct(LAND_TYPE), &
          & parent=box_tile, &
          & fract_varname='equal')   ! Fraction = 1.
      END IF

        IF (model%config%use_glacier) THEN
          glacier_tile => t_jsb_tile('glac', 'Glacier tile', processes=(/SSE_/), &
            & lct=t_jsb_lct(GLACIER_TYPE), parent=land_tile, &
            & fract_varname='fract_glac')
        END IF

          veg_tile => t_jsb_tile('veg', 'Veg tile', processes=(/SSE_, PHENO_/), &
            & lct=t_jsb_lct(VEG_TYPE), parent=land_tile,  &
            & fract_varname='fract_veg')

    NULLIFY(box_tile, land_tile, veg_tile)
    IF (model%config%use_glacier) NULLIFY(glacier_tile)
    IF (model%config%use_lakes) NULLIFY(lake_tile)

  END SUBROUTINE init_usecase_lite

  SUBROUTINE init_usecase_lite_new(model)

    TYPE(t_jsb_model), POINTER, INTENT(inout) :: model

    CLASS(t_jsb_tile), POINTER :: box_tile, lake_tile, land_tile, glacier_tile, veg_tile

    CHARACTER(len=*), PARAMETER :: routine = modname//':init_usecase_lite_new'

    ! Create root of the tile structure representing the whole grid box
    box_tile => t_jsb_tile('box', 'Top tile', &
      ! & processes      =(/A2L_,     L2A_,     RAD_,      HYDRO_,    TURB_    /), &
      ! & process_actions=(/ON_TILE_, ON_TILE_, ON_LEAFS_, ON_LEAFS_, ON_LEAFS_/), &
      & processes      =(/A2L_,     L2A_,     SEB_,      RAD_,      HYDRO_,    TURB_    /), &
      & process_actions=(/ON_TILE_, ON_TILE_, ON_LEAFS_, ON_LEAFS_, ON_LEAFS_, ON_LEAFS_/), &
      & model_id=model%id, fract_varname='notsea')
    model%top => box_tile

#ifndef __NO_JSBACH_HD__
    IF (model%processes(HD_)%p%config%active) THEN
      CALL box_tile%Set_process_action(HD_, ON_TILE_)
    END IF
#endif

      ! Create sub-tiles
      IF (model%config%use_lakes) THEN
        lake_tile => t_jsb_tile('lake', 'Lake tile', &
          ! & processes=(/SEB_/), process_actions=(/ON_TILE_/), &
          & lct=t_jsb_lct(LAKE_TYPE), parent=box_tile, fract_varname='fract_lake')
        land_tile => t_jsb_tile('land', 'Land tile', &
          ! & processes=(/SEB_/), process_actions=(/ON_TILE_/), &
          & processes=(/SSE_/), process_actions=(/ON_LEAFS_/), &
          & lct=t_jsb_lct(LAND_TYPE), parent=box_tile)
      ELSE
        land_tile => t_jsb_tile('land', 'Land tile', &
          ! & processes=(/SEB_/), process_actions=(/ON_TILE_/), &
          & processes=(/SSE_/), process_actions=(/ON_LEAFS_/), &
          & lct=t_jsb_lct(LAND_TYPE), &
          & parent=box_tile, &
          & fract_varname='equal')   ! Fraction = 1.
      END IF

        IF (model%config%use_glacier) THEN
          glacier_tile => t_jsb_tile('glac', 'Glacier tile', &
            & lct=t_jsb_lct(GLACIER_TYPE), parent=land_tile, &
            & fract_varname='fract_glac')

          veg_tile => t_jsb_tile('veg', 'Veg tile', &
            & processes=(/PHENO_/), process_actions=(/ON_TILE_/), &
            & lct=t_jsb_lct(VEG_TYPE), parent=land_tile,  &
            & fract_varname='fract_veg')
        ELSE
          veg_tile => t_jsb_tile('veg', 'Veg tile', &
            & processes=(/PHENO_/), process_actions=(/ON_TILE_/), &
            & lct=t_jsb_lct(VEG_TYPE), parent=land_tile,  &
            & fract_varname='equal')
        END IF

    NULLIFY(box_tile, land_tile, veg_tile)
    IF (model%config%use_glacier) NULLIFY(glacier_tile)
    IF (model%config%use_lakes) NULLIFY(lake_tile)

  END SUBROUTINE init_usecase_lite_new

  SUBROUTINE init_usecase_pfts(model)

    TYPE(t_jsb_model), POINTER, INTENT(inout) :: model

    CLASS(t_jsb_tile), POINTER :: box_tile, lake_tile, land_tile, glacier_tile, veg_tile, pft_tile

    INTEGER :: i

    INTEGER, PARAMETER :: npft = 11

    CHARACTER(len=5), PARAMETER :: pft_shortnames(npft) = [character(len=5) ::  &
      & 'pft01', 'pft02', 'pft03', 'pft04', 'pft05', 'pft06', &
      & 'pft07', 'pft08', 'pft09', 'pft10', 'pft11' ]

    CHARACTER(len=30), PARAMETER :: pft_longnames(npft) = [character(len=30) :: &
      & 'Tropical broadleaf evergreen', &
      & 'Tropical broadleaf deciduous', &
      & 'Extra-tropical evergreen',     &
      & 'Extra-tropical deciduous',     &
      & 'Raingreen shrubs',             &
      & 'Deciduous shrubs',             &
      & 'C3 grass',                     &
      & 'C4 grass',                     &
      & 'C3 pasture',                   &
      & 'C4 pasture',                   &
      & 'C3 crops'                    ]

    ! Index of PFT into lctlib columns
    ! id is one less than in lctlib file because the first column (glacier) has been thrown away
    INTEGER, PARAMETER :: pft_ids(npft) = (/ 1,2,3,4,9,10,11,12,14,15,19 /)

    CHARACTER(len=*), PARAMETER :: routine = modname//':init_usecase_pfts'


    ! Create root of the tile structure representing the whole grid box
    box_tile => t_jsb_tile('box', 'Top tile', &
      & processes      =(/A2L_,     L2A_,     RAD_  ,    TURB_    /), &
      & process_actions=(/ON_TILE_, ON_TILE_, ON_LEAFS_, ON_LEAFS_/), &
      & model_id=model%id, fract_varname='notsea')
    model%top => box_tile

#ifndef __NO_JSBACH_HD__
    IF (model%processes(HD_)%p%config%active) THEN
      CALL box_tile%Set_process_action(HD_, ON_TILE_)
    END IF
#endif

    IF (model%processes(PPLCC_)%p%config%active) THEN
      CALL box_tile%Set_process_action(PPLCC_, ON_LEAFS_)
    END IF

      IF (model%config%use_lakes) THEN
        lake_tile => t_jsb_tile('lake', 'Lake tile', &
          & processes=(/SEB_, HYDRO_/), process_actions=(/ON_TILE_, ON_TILE_/), &
          & lct=t_jsb_lct(LAKE_TYPE), parent=box_tile)
        land_tile => t_jsb_tile('land', 'Land tile', &
          & processes=(/SEB_/), process_actions=(/ON_TILE_/), &
          & lct=t_jsb_lct(LAND_TYPE), parent=box_tile)
      ELSE
        land_tile => t_jsb_tile('land', 'Land tile', &
          & processes=(/SEB_/), process_actions=(/ON_TILE_/), &
          & lct=t_jsb_lct(LAND_TYPE), parent=box_tile, &
          & fract_varname='equal')   ! Fraction = 1.
      END IF

        IF (model%config%use_glacier) THEN
          glacier_tile => t_jsb_tile('glac', 'Glacier tile', &
            & processes = (/HYDRO_, SSE_/), &
            & process_actions = (/ON_TILE_, ON_TILE_/), &
            & lct=t_jsb_lct(GLACIER_TYPE), parent=land_tile)
            ! & lct=t_jsb_lct(GLACIER_TYPE, lib_id=1), parent=land_tile)
        END IF

        ! default jsb4: HYDRO_ is active at veg tile only to save runtime
        IF (.NOT. model%config%use_quincy) THEN
          veg_tile => t_jsb_tile('veg', 'Veg tile', &
            & processes       = (/HYDRO_,   SSE_,     ASSIMI_,   PHENO_/),    &
            & process_actions = (/ON_TILE_, ON_TILE_, ON_LEAFS_, ON_LEAFS_/), &
            & lct=t_jsb_lct(VEG_TYPE), parent=land_tile,  &
            & fract_varname='fract_veg')

          IF (model%processes(ALCC_)%p%config%active) THEN
            CALL veg_tile%Set_process_action(ALCC_, ON_TILE_, AGGREGATE_)
          END IF
          !
          IF (model%processes(NLCC_)%p%config%active) THEN
            CALL veg_tile%Set_process_action(NLCC_, ON_TILE_, AGGREGATE_)
          END IF
          !
          ! use_quincy: 
        !  HYDRO_ is active at PFT tiles to make all the HYDRO_ variables available to the PFT tiles
        !    note: HYDRO_ results differ only marginally between ON_TILE_ & ON_LEAFS_ (tested by Reiner)
        !  VEG_ is active
        ELSE
          veg_tile => t_jsb_tile('veg', 'Veg tile', &
            & processes       = (/HYDRO_,    SSE_,     ASSIMI_,   PHENO_,    VEG_   /), &
            & process_actions = (/ON_LEAFS_, ON_TILE_, ON_LEAFS_, ON_LEAFS_, ON_LEAFS_/), &
            & lct=t_jsb_lct(VEG_TYPE), parent=land_tile,  &
            & fract_varname='fract_veg')
        ENDIF

        IF (model%processes(TLCC_)%p%config%active) THEN
          CALL veg_tile%Set_process_action(TLCC_, ON_TILE_, AGGREGATE_)
        END IF

          ! R: TBD: DISTURB_ ist kein Task von CARBON_,
          !         da es die cover fractions der pfts verändert und nur als Nebeneffekt auch die C-Pools
          !         verändert. Er würde im Prinzip auch ablaufen ohne den Prozess CARBON_ und würde auch
          !         auf alle anderen Stoffpools wirken. Außerdem soll CARBON_ später wie NITROGEN_
          !         oder andere Stoffpools gehandelt werden. Diese Prozesse haben bis
          !         auf ihre Fixierung kaum eigene Tasks sondern werden von anderen Prozessen verändert.
          !         => die Stoffpools sollten langfristig evtl. gar keine Prozesse sein
          !
          !         DISTURB_ braucht CARBON,
          !         da es C-Pools als Input für die Feuerberechnung braucht und sollte diese auch verändern
          !         wenn der CARBON Prozess existiert. Dies ist ein grundsätzliches Problem: Prozesse benötigen
          !         - oder sogar verändern - evtl. Variablen, die eigentlich zu einem anderen Prozess gehören.
          !         Das bedeutet der Speicher des anderen Prozesses muß angelegt werden auch wenn dieser
          !         Prozess gar nicht abläuft. Dazu müsste wohl am besten eine Prozedur die die Namelist-
          !         Einträge entwprechend umändert geschrieben werden oder wenn Usecases benutzt werden
          !         (wie hier) der entsprechende Prozess dennoch definiert werden, auch wenn in der Queue-
          !         Liste nachher kein Task dieses Prozesses auftaucht.
          !
          !         DISTURB_ braucht PHENO wg. veg_fract_correct s.o.
          !
          !         CARBON_ braucht DISTURB_ nicht!
          !
          !         Manche Prozeduren sollen nur ausgeführt werden wenn beide Prozesse aktiv sind.
          !         In DISTURB_ z.B. die dazugehörige C-Umverteilung bei Feuer und Winbruch.
          !         Dies habe ich mit "IF ( tile%Is_process_active(CARBON_) ) THEN" gemacht.
          !
          IF (model%processes(CARBON_)%p%config%active) THEN
            CALL veg_tile%Set_process_action(CARBON_, ON_LEAFS_, AGGREGATE_)
          END IF
          IF (model%processes(TCQ_)%p%config%active) THEN
            CALL veg_tile%Set_process_action(TCQ_, ON_LEAFS_, AGGREGATE_)
          END IF
          IF (model%processes(FUEL_)%p%config%active) THEN
            CALL veg_tile%Set_process_action(FUEL_, ON_TILE_, AGGREGATE_)
          END IF
          IF (model%processes(DISTURB_)%p%config%active) THEN
            CALL veg_tile%Set_process_action(DISTURB_, ON_LEAFS_, AGGREGATE_)
          END IF

            DO i=1,npft
              pft_tile => t_jsb_tile(pft_shortnames(i), TRIM(pft_longnames(i)), &
               & lct=t_jsb_lct(VEG_TYPE, lib_id=pft_ids(i)), parent=veg_tile)
            END DO

    NULLIFY(box_tile, land_tile, veg_tile, pft_tile)
    IF (model%config%use_glacier) NULLIFY(glacier_tile)
    IF (model%config%use_lakes) NULLIFY(lake_tile)

  END SUBROUTINE init_usecase_pfts

#endif
END MODULE mo_jsb_model_usecases
