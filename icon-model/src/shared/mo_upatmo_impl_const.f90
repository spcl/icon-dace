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

! Constants comparable to those in 'src/shared/mo_impl_constants'.

MODULE mo_upatmo_impl_const

  USE mo_kind,    ONLY: wp

  IMPLICIT NONE

  PRIVATE

  ! General:
  PUBLIC :: imsg_thr
  PUBLIC :: itmr_thr
  PUBLIC :: iUpatmoStat
  ! Physics:
  PUBLIC :: isolvar
  PUBLIC :: isolvardat
  PUBLIC :: iorbit
  PUBLIC :: icycle
  PUBLIC :: iUpatmoPrcId
  PUBLIC :: iUpatmoGrpId
  PUBLIC :: iUpatmoPrcMode
  PUBLIC :: iUpatmoTendId
  PUBLIC :: iUpatmoPrcStat
  PUBLIC :: iUpatmoGasId
  PUBLIC :: iUpatmoGasMode
  PUBLIC :: iUpatmoGasStat
  PUBLIC :: iUpatmoExtdatId
  PUBLIC :: iUpatmoExtdatStat
  PUBLIC :: iUpatmoExtdatLatId
  PUBLIC :: iUpatmoExtdatLevId
  PUBLIC :: iUpatmoExtdatTimeId
  PUBLIC :: iUpatmoTracerId
  PUBLIC :: startHeightDef
  PUBLIC :: iThermdynCoupling

  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------
  !                                    Notes
  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

  ! Some of the constants below (indicated by "(nitem)") 
  ! follow the construction rule:
  !
  ! TYPE t_ixyz
  !   INTEGER :: a      ! Identifier 1
  !   INTEGER :: b      ! Identifier 2
  !   INTEGER :: c      ! Identifier 3
  !   INTEGER :: d      ! Identifier 4
  !   ! ...
  !   !
  !   INTEGER :: nitem  ! Number of identifiers in t_ixyz
  ! END TYPE t_ixyz
  ! TYPE(t_ixyz), PARAMETER :: ixyz = t_ixyz( 1, & ! a
  !                                           2, & ! b = a + 1
  !                                           3, & ! c = b + 1
  !                                           4, & ! d = c + 1
  !
  !                                           4  ) ! nitem = d, since d is the last list element
  !
  ! Possible advantages of this construction rule:
  ! * We may want to loop over whatever the items of ixyz identify: 
  !
  !    DO jitem = 1, ixyz%nitem
  !      ! jitem successivly takes the values of ixyz%a, ixyz%b, ixyz%c, ... 
  !      ! and can be compared with namelist entries etc.
  !    ENDDO
  ! 
  ! Disadvantages of this construction rule:
  ! * The identifiers, selectable by namelist are fixed to: 1, 2, 3, ... 
  !   (something like: 5, 6, 8, 20, 124, ... is not possible).  
  !   An array xyz(1:ixyz%nitem), which contains the actual identifiers, 
  !   selectable by namelist could handle this, but adds a further level 
  !   of complexity, which we would like to avoid.
  !
  ! Please, be careful, if you modify the following types, 
  ! since a reasonable check for the adherence to the above construction rule 
  ! during runtime is not possible.

  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------
  !                                  General
  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

  !-----------------------------------------------------
  !                    Thresholds
  !-----------------------------------------------------

  TYPE t_ithr
    INTEGER :: low    ! Low threshold
    INTEGER :: med    ! Medium threshold
    INTEGER :: high   ! High threshold
  END TYPE t_ithr
  !
  ! ... for message output:
  TYPE(t_ithr), PARAMETER :: imsg_thr = t_ithr(  8, &  !imsg_thr%low
    &                                           10, &  !imsg_thr%med
    &                                           15  )  !imsg_thr%high
  !
  ! ... for timers:
  TYPE(t_ithr), PARAMETER :: itmr_thr = t_ithr(  2, &  !itmr_thr%low
    &                                            5, &  !itmr_thr%med
    &                                            8  )  !itmr_thr%high

  !-----------------------------------------------------
  !             Configuration status (nitem)
  !-----------------------------------------------------

  TYPE t_iUpatmoStat
    INTEGER :: checked      ! Upper-atmosphere namelist settings crosschecked?
    INTEGER :: configured   ! Upper-atmosphere configured?
    INTEGER :: required     ! Upper-atmosphere settings required at all?
    INTEGER :: message      ! Message output desired?
    INTEGER :: timer        ! Timer monitoring desired?
    !
    INTEGER :: nitem        ! Number of identifiers
  END TYPE t_iUpatmoStat
  TYPE(t_iUpatmoStat), PARAMETER :: iUpatmoStat = t_iUpatmoStat( 1, &  !iUpatmoStat%checked
    &                                                            2, &  !iUpatmoStat%configured
    &                                                            3, &  !iUpatmoStat%required
    &                                                            4, &  !iUpatmoStat%message
    &                                                            5, &  !iUpatmoStat%timer
    !
    &                                                            5  )  !iUpatmoStat%nitem

  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------
  !                           Upper-atmosphere physics
  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

  !-----------------------------------------------------
  !               Solar activity (nitem)
  !-----------------------------------------------------

  TYPE t_isolvar
    INTEGER :: norm    ! Normal activity
    INTEGER :: low     ! Low --,,--
    INTEGER :: high    ! High --,,--
    !
    INTEGER :: nitem   ! Number of entries
  END TYPE t_isolvar
  TYPE(t_isolvar), PARAMETER :: isolvar = t_isolvar( 1, &  !isolvar%norm
    &                                                2, &  !isolvar%low
    &                                                3, &  !isolvar%high
    !
    &                                                3  )  !isolvar%nitem

  !-----------------------------------------------------
  !             Solar activity data (nitem)
  !-----------------------------------------------------

  TYPE t_isolvardat
    INTEGER :: rottman   ! G. Rottman data
    INTEGER :: lean      ! J. Lean data
    !
    INTEGER :: nitem     ! Number of entries
  END TYPE t_isolvardat
  TYPE(t_isolvardat), PARAMETER :: isolvardat = t_isolvardat( 1, &  !isolvardat%rottman
    &                                                         2, &  !isolvardat%lean
    !
    &                                                         2  )  !isolvardat%nitem

  !-----------------------------------------------------
  !                Orbit model (nitem)
  !-----------------------------------------------------

  TYPE t_iorbit
    INTEGER :: vsop87  ! Standard and accurate model
    INTEGER :: kepler  ! Simple model, appropriate for idealized work 
    !
    INTEGER :: nitem   ! Number of entries
  END TYPE t_iorbit
  TYPE(t_iorbit), PARAMETER :: iorbit = t_iorbit( 1, &  !iorbit%vsop87
    &                                             2, &  !iorbit%kepler
    !
    &                                             2  )  !iorbit%nitem

  !-----------------------------------------------------
  !                Solar cycle (nitem)
  !-----------------------------------------------------

  TYPE t_icycle
    INTEGER :: std     ! Standard setting
    INTEGER :: day27   ! 27day cycle
    !
    INTEGER :: nitem   ! Number of entries
  END TYPE t_icycle
  TYPE(t_icycle), PARAMETER :: icycle = t_icycle( 1, &  !icycle%std
    &                                             2, &  !icycle%day27
    !
    &                                             2  )  !icycle%nitem

  !------------------------------------------------------
  !      Identifiers for physics processes (nitem)
  !------------------------------------------------------

  TYPE t_iUpatmoPrcId
    INTEGER :: vdfmol   ! Molecular diffusion
    INTEGER :: fric     ! Frictional heating
    INTEGER :: iondrag  ! Ion drag
    INTEGER :: joule    ! Joule heating
    INTEGER :: srbc     ! Heating due to Schumann-Runge bands and continuum of O2
    INTEGER :: nlte     ! Heating due to processes in non-localc-thermal-equilibrium
    INTEGER :: euv      ! Extreme-ultraviolet heating
    INTEGER :: no       ! Near-infrared heating by NO
    INTEGER :: chemheat ! Chemical heating
    !
    INTEGER :: nitem    ! Number of entries
  END TYPE t_iUpatmoPrcId
  TYPE(t_iUpatmoPrcId), PARAMETER :: iUpatmoPrcId = t_iUpatmoPrcId( 1, &  !iUpatmoPrcId%vdfmol
    &                                                               2, &  !iUpatmoPrcId%fric
    &                                                               3, &  !iUpatmoPrcId%iondrag
    &                                                               4, &  !iUpatmoPrcId%joule
    &                                                               5, &  !iUpatmoPrcId%srbc
    &                                                               6, &  !iUpatmoPrcId%nlte
    &                                                               7, &  !iUpatmoPrcId%euv
    &                                                               8, &  !iUpatmoPrcId%no
    &                                                               9, &  !iUpatmoPrcId%chemheat
    !
    &                                                               9  )  !iUpatmoPrcId%nitem

  !------------------------------------------------------
  !     Default start heights for physics processes
  !------------------------------------------------------

  TYPE t_startHeightDef
    REAL(wp) :: vdfmol    ! For molecular diffusion
    REAL(wp) :: fric      ! For frictional heating
    REAL(wp) :: iondrag   ! For ion drag
    REAL(wp) :: joule     ! For joule heating
    REAL(wp) :: srbc      ! For heating due to Schumann-Runge bands and continuum of O2
    REAL(wp) :: nlte      ! For heating due to processes in non-localc-thermal-equilibrium
    REAL(wp) :: euv       ! For extreme-ultraviolet heating
    REAL(wp) :: no        ! For near-infrared heating by NO
    REAL(wp) :: chemheat  ! For chemical heating
  END TYPE t_startHeightDef
  TYPE(t_startHeightDef), PARAMETER :: startHeightDef = t_startHeightDef( 75000._wp, &  ! (m) startHeightDef%vdfmol
    &                                                                     85000._wp, &  ! (m) startHeightDef%fric
    &                                                                     80000._wp, &  ! (m) startHeightDef%iondrag
    &                                                                     80000._wp, &  ! (m) startHeightDef%joule = 
    !                                                                                         startHeightDef%iondrag !
    &                                                                     50000._wp, &  ! (m) startHeightDef%srbc
    &                                                                         0._wp, &  ! (m) startHeightDef%nlte
    &                                                                     90000._wp, &  ! (m) startHeightDef%euv
    &                                                                     60000._wp, &  ! (m) startHeightDef%no
    &                                                                     70000._wp  )  ! (m) startHeightDef%chemheat
  
  !-----------------------------------------------------
  !    Identifiers of physics process groups (nitem)
  !-----------------------------------------------------

  TYPE t_iUpatmoGrpId
    INTEGER :: imf     ! Ion drag (I), molecular diffusion (M) and frictional heating (F)
    INTEGER :: rad     ! Radiation and chemical heating 
    !
    INTEGER :: nitem   ! Number of entries
  END TYPE t_iUpatmoGrpId
  TYPE(t_iUpatmoGrpId), PARAMETER :: iUpatmoGrpId = t_iUpatmoGrpId( 1, &  !iUpatmoGrpId%imf
    &                                                               2, &  !iUpatmoGrpId%rad
    !
    &                                                               2  )  !iUpatmoGrpId%nitem

  !-----------------------------------------------------
  !     Identifiers for modes of physics process 
  !                and physics groups
  !-----------------------------------------------------

  TYPE t_iUpatmoPrcMode
    INTEGER :: unassigned  ! For configuration purposes
    INTEGER :: off         ! Process is disabled
    INTEGER :: on          ! Process is enabled
    INTEGER :: offline     ! Tendencies are computed, but not applied
    !
    INTEGER :: startitem   ! Integer of firt list entry (since it differes from 1)
    INTEGER :: enditem     ! Integer of last list entry 
                           ! (not nitem, since it differs from number of entries)
  END TYPE t_iUpatmoPrcMode
  TYPE(t_iUpatmoPrcMode), PARAMETER :: iUpatmoPrcMode = t_iUpatmoPrcMode( -1, &  !iUpatmoPrcMode%unassigned
    &                                                                      0, &  !iUpatmoPrcMode%off
    &                                                                      1, &  !iUpatmoPrcMode%on
    &                                                                      2, &  !iUpatmoPrcMode%offline
    !
    &                                                                     -1, &  !iUpatmoPrcMode%startitem
    &                                                                      2  )  !iUpatmoPrcMode%enditem

  !-----------------------------------------------------
  !        Identifiers for variables for which 
  !    physics processes compute tendencies (nitem)
  !-----------------------------------------------------

  TYPE t_iUpatmoTendId
    INTEGER :: temp     ! Temperature
    INTEGER :: wind_h   ! Zonal and meridional wind components u and v
    INTEGER :: qx       ! Tracer (water vapor)
    !
    INTEGER :: nitem    ! Number of entries up to here
    ! Derived
    INTEGER :: exner    ! Exner pressure
    !
    INTEGER :: nitem_2  ! Number of entries up to here without 'nitem'
  END TYPE t_iUpatmoTendId
  TYPE(t_iUpatmoTendId), PARAMETER :: iUpatmoTendId = t_iUpatmoTendId( 1, &  !iUpatmoTendId%temp
    &                                                                  2, &  !iUpatmoTendId%wind_h
    &                                                                  3, &  !iUpatmoTendId%qx
    !
    &                                                                  3, &  !iUpatmoTendId%nitem
    !
    &                                                                  4, &  !iUpatmoTendId%exner
    !
    &                                                                  4  )  !iUpatmoTendId%nitem_2


  !-----------------------------------------------------
  !           Identifiers for physics status 
  !           (in the broadest sense) (nitem)
  !-----------------------------------------------------

  TYPE t_iUpatmoPrcStat
    INTEGER :: enabled          ! Parameterization is switched on
    INTEGER :: offline          ! Compute, but not apply tendencies
    INTEGER :: initialized      ! Parameterization is initialized
    INTEGER :: finalized        ! Parameterization is finalized
    INTEGER :: afterActivePhase ! Active phase of parameterization is over
    ! 
    INTEGER :: nitem        ! Number of entries
  END TYPE t_iUpatmoPrcStat
  TYPE(t_iUpatmoPrcStat), PARAMETER :: iUpatmoPrcStat = t_iUpatmoPrcStat( 1, &  !iUpatmoPrcStat%enabled
    &                                                                     2, &  !iUpatmoPrcStat%offline
    &                                                                     3, &  !iUpatmoPrcStat%initialized
    &                                                                     4, &  !iUpatmoPrcStat%finalized
    &                                                                     5, &  !iUpatmoPrcStat%afterActivePhase
    !
    &                                                                     5  )  !iUpatmoPrcStat%nitem

  !-----------------------------------------------------
  !   Identifiers of radiatively active gases (nitem)
  !-----------------------------------------------------

  TYPE t_iUpatmoGasId
    INTEGER :: o3      ! Ozone
    INTEGER :: o2      ! Dioxygen
    INTEGER :: o       ! Atomar oxygen
    INTEGER :: co2     ! Carbon dioxide
    INTEGER :: no      ! Nitric oxide
    INTEGER :: n2      ! Dinitrogen
    ! 
    INTEGER :: nitem   ! Number of entries
    ! Selectors
    INTEGER :: diag    ! Which gas concentration might be computed
                       ! diagnostically from all other gas concentrations    
  END TYPE t_iUpatmoGasId
  TYPE(t_iUpatmoGasId), PARAMETER :: iUpatmoGasId = t_iUpatmoGasId( 1, &  !iUpatmoGasId%o3
    &                                                               2, &  !iUpatmoGasId%o2
    &                                                               3, &  !iUpatmoGasId%o
    &                                                               4, &  !iUpatmoGasId%co2
    &                                                               5, &  !iUpatmoGasId%no
    &                                                               6, &  !iUpatmoGasId%n2
    !
    &                                                               6, &  !iUpatmoGasId%nitem
    !
    &                                                               6  )  !iUpatmoGasId%diag -> n2

  !-----------------------------------------------------
  !             Identifiers for gas mode
  !-----------------------------------------------------

  TYPE t_iUpatmoGasMode
    INTEGER :: zero        ! Zero gas concentration
    INTEGER :: const       ! Horizontally/vertically/temporally constant 
                           ! (single fixed value, read from namelist)          
    INTEGER :: extdat      ! External data read from file
    INTEGER :: diag        ! (Only for N2): determine N2 as the residual of all other gases
    ! 
    INTEGER :: startitem   ! Integer of firt list entry (since it differes from 1)
    INTEGER :: enditem     ! Integer of last list entry 
                           ! (not nitem, since it differs from number of entries)    
  END TYPE t_iUpatmoGasMode
  TYPE(t_iUpatmoGasMode), PARAMETER :: iUpatmoGasMode = t_iUpatmoGasMode( 0, &  !iUpatmoGasMode%zero
    &                                                                     1, &  !iUpatmoGasMode%const
    &                                                                     2, &  !iUpatmoGasMode%extdat
    &                                                                     3, &  !iUpatmoGasMode%diag
    !
    &                                                                     0, &  !iUpatmoGasMode%startitem
    &                                                                     3  )  !iUpatmoGasMode%enditem
  
  !-----------------------------------------------------
  !         Identifiers for gas status (nitem)
  !-----------------------------------------------------

  TYPE t_iUpatmoGasStat
    INTEGER :: enabled      ! Gas is required
    INTEGER :: initialized  ! Gas is initialized
    INTEGER :: finalized    ! Gas is finalized
    !
    INTEGER :: nitem        ! Number of entries
  END TYPE t_iUpatmoGasStat
  TYPE(t_iUpatmoGasStat), PARAMETER :: iUpatmoGasStat = t_iUpatmoGasStat( 1, &  !iUpatmoGasStat%enabled
    &                                                                     2, &  !iUpatmoGasStat%initialized
    &                                                                     3, &  !iUpatmoGasStat%finalized
    !
    &                                                                     3  )  !iUpatmoGasStat%nitem

  !-----------------------------------------------------
  !       Identifiers of external data (nitem)
  !-----------------------------------------------------

  TYPE t_iUpatmoExtdatId
    INTEGER :: gases      ! Radiatively active gases
    INTEGER :: chemheat   ! Chemical heating
    ! 
    INTEGER :: nitem      ! Number of entries
  END TYPE t_iUpatmoExtdatId
  TYPE(t_iUpatmoExtdatId), PARAMETER :: iUpatmoExtdatId = t_iUpatmoExtdatId( 1, &  !iUpatmoExtdatId%gases
    &                                                                        2, &  !iUpatmoExtdatId%chemheat
    !
    &                                                                        2  )  !iUpatmoExtdatId%nitem

  !-----------------------------------------------------
  !        Identifiers of extdat status (nitem)
  !-----------------------------------------------------

  TYPE t_iUpatmoExtdatStat
    INTEGER :: required     ! External data is required
    INTEGER :: initialized  ! External data is initialized
    INTEGER :: finalized    ! External data is finalized
    ! 
    INTEGER :: nitem      ! Number of entries
  END TYPE t_iUpatmoExtdatStat
  TYPE(t_iUpatmoExtdatStat), PARAMETER :: iUpatmoExtdatStat = t_iUpatmoExtdatStat( 1, &  !iUpatmoExtdatStat%required
    &                                                                              2, &  !iUpatmoExtdatStat%initialized
    &                                                                              3, &  !iUpatmoExtdatStat%finalized
    !
    &                                                                              3  )  !iUpatmoExtdatStat%nitem

  !------------------------------------------------------
  !   Identifiers of external data: latitude (nitem)
  !------------------------------------------------------

  TYPE t_iUpatmoExtdatLatId
    INTEGER :: deg        ! Latitudes in degree north
    ! 
    INTEGER :: nitem      ! Number of entries
  END TYPE t_iUpatmoExtdatLatId
  TYPE(t_iUpatmoExtdatLatId), PARAMETER :: iUpatmoExtdatLatId = t_iUpatmoExtdatLatId( 1, &  !iUpatmoExtdatLatId%deg
    !
    &                                                                                 1  )  !iUpatmoExtdatLatId%nitem

  !-----------------------------------------------------
  !    Identifiers of external data: level (nitem)
  !-----------------------------------------------------

  TYPE t_iUpatmoExtdatLevId
    INTEGER :: p          ! Pressure levels
    INTEGER :: z          ! Geometric height levels
    ! 
    INTEGER :: nitem      ! Number of entries
  END TYPE t_iUpatmoExtdatLevId
  TYPE(t_iUpatmoExtdatLevId), PARAMETER :: iUpatmoExtdatLevId = t_iUpatmoExtdatLevId( 1, &  !iUpatmoExtdatLevId%p
    &                                                                                 2, &  !iUpatmoExtdatLevId%z
    !
    &                                                                                 2  )  !iUpatmoExtdatLevId%nitem

  !------------------------------------------------------
  !     Identifiers of external data: time (nitem)
  !------------------------------------------------------

  TYPE t_iUpatmoExtdatTimeId
    INTEGER :: month      ! Times in months
    ! 
    INTEGER :: nitem      ! Number of entries
  END TYPE t_iUpatmoExtdatTimeId
  TYPE(t_iUpatmoExtdatTimeId), PARAMETER :: iUpatmoExtdatTimeId = t_iUpatmoExtdatTimeId( 1, &  !iUpatmoExtdatTimeId%month
    !
    &                                                                                    1  )  !iUpatmoExtdatTimeId%nitem

  !----------------------------------------------------------
  ! Identifiers of upper-atmosphere-affected tracers (nitem)
  !----------------------------------------------------------

  ! Please note that they are not (necessarily) identical 
  ! to the 'iqx' in 'src/configure_model/mo_run_config', 
  ! set in 'src/configure_model/mo_nml_crosscheck'.

  TYPE t_iUpatmoTracerId
    INTEGER :: qv         ! Water vapor
    ! 
    INTEGER :: nitem      ! Number of entries
  END TYPE t_iUpatmoTracerId
  TYPE(t_iUpatmoTracerId), PARAMETER :: iUpatmoTracerId = t_iUpatmoTracerId( 1, &  !iUpatmoTracerId%qv
    !
    &                                                                        1  )  !iUpatmoTracerId%nitem  

  !----------------------------------------------------------
  !      Identifiers for thermodynamic coupling between
  !               physics and dynamics (nitem)
  !----------------------------------------------------------

  TYPE t_iThermdynCoupling
    INTEGER :: isobaric   ! coupling assumes constant pressure
    INTEGER :: isochoric  ! coupling assumes constant volume/mass
    INTEGER :: entropic   ! coupling via heat (entropy) source
    ! 
    INTEGER :: nitem      ! Number of entries
  END TYPE t_iThermdynCoupling
  TYPE(t_iThermdynCoupling), PARAMETER :: iThermdynCoupling = t_iThermdynCoupling( 1, &  !iThermdynCoupling%isobaric
    &                                                                              2, &  !iThermdynCoupling%isochoric
    &                                                                              3, &  !iThermdynCoupling%entropic
    !
    &                                                                              3  )  !iThermdynCoupling%nitem

END MODULE mo_upatmo_impl_const
