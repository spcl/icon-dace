! mo_jsbach_to_hd.f90 - Various routines for the subroutine couping of JSBACH to HD
!
! Copyright (C) 2014, MPI-M
! SPDX-License-Identifier: BSD-3-Clause
! See ./LICENSES/ for license information
!_________________________________________

MODULE mo_jsbach_to_hd
  !
  !
  ! Note that the coupled 0.5 degree version is also part of MPI-ESM, the Earth System Model of the 
  !    Max Planck Institute for Meteorology (Mauritsen et al. 2019). 
  !    Reference: Mauritsen, T., et al. (2019) Developments in the MPI-M Earth System Model 
  !    version 1.2 (MPI-ESM1.2) and its response to increasing CO2. J. Adv. Model. Earth Syst., 11, 
  !    doi: 10.1029/2018MS001400.
  !
  ! This module comprises several routines that are used within MPI-ESM where HD model is called
  ! as a subroutine of the land surface model JSBACH. As these are not used, when HD is run as
  ! an independent model, they were taken out from their original modules. A usage within MPI-ESM 
  ! has not been tested (July 2022). 
  ! In addition, it currently has not been included in the compile and dependency files.
  !
  ! Source            Subroutine        
  ! ------            ----------        
  ! mo_hydrology
  !    Public:        hydrology_collect 
  !    Others:        hydrology_corr, hydrology_echam
  !                   intpol, intpol_compute_mapping, intpol_coord_axis_setup
  !                   intpol_with_mapping, prepare_intpol_mapping, reassign_runoff
  !    New routines with code parts:
  !                   init_jsbach_to_hd, cleanup_jsbach_to_hd
  !                   
  ! 
  !

  USE mo_array_utils,      ONLY: dec_monotonic_closest_midpoint
  USE mo_constants,        ONLY: rhoh2o
  USE mo_control,          ONLY: nlon, ngl
  USE mo_gaussgrid,        ONLY: philat, philon
  USE mo_grid,             ONLY: domain
  USE mo_hydrology,        ONLY: grid_hd, hd_calling_interval
  USE mo_kind,             ONLY: dp
  USE mo_mpi,              ONLY: p_parallel_io 
  USE mo_time_control,     ONLY: l_puthd, l_gethd, puthd, gethd, delta_time,  &
                                 ev_puthd, get_interval_seconds_next

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: hydrology_collect

  TYPE cart_idx_2d
    INTEGER :: ilon, ilat
  END TYPE cart_idx_2d

  TYPE(cart_idx_2d),  ALLOCATABLE :: intpol_mapping(:,:) ! no longer used (replaced by scrip remapping)
  REAL(dp), PARAMETER :: fullcirc = 360.0_dp



CONTAINS
! 
! ********** New Routines ***************************************

  SUBROUTINE init_jsbach_to_hd(nl,nb)

    INTEGER,    INTENT(in)    :: nl      ! Number of HD longitudes
    INTEGER,    INTENT(in)    :: nb      ! Number of HD latitudes

    IF (p_parallel_io) THEN
      ALLOCATE (intpol_mapping(nl, nb)) ; intpol_mapping = cart_idx_2d(-1, -1)
      hd_calling_interval = get_interval_seconds_next(ev_puthd)
    ENDIF

  END SUBROUTINE init_jsbach_to_hd

  SUBROUTINE cleanup_jsbach_to_hd

    IF (p_parallel_io) THEN
      DEALLOCATE  (intpol_mapping)
    ENDIF

  END SUBROUTINE cleanup_jsbach_to_hd
! 
! ********** Coupling Routines ***************************************

  SUBROUTINE hydrology_collect (knproma,            &
                                paros, padrain,     &
                                papmecal,           &
                                pdisch, pruntoc,    &
                                pros_hd, pdrain_hd, &
                                palac)

    !  Collects runoff and drainage  for input to the HD-model
    !
    !  hydrology_collect is called from physc
    !
    ! Authors:
    !
    ! U. Schlese, MPI, August 2000
    ! I. Kirchner, MPI, April 2001, date/time control
    ! L. Kornblueh, MPI, July 2002, parallelization

    !  scalar arguments

    INTEGER, INTENT(in) :: knproma

    ! array arguments

    REAL(dp), INTENT(inout) :: paros(knproma),   &! acc. runoff for HD-model
                               padrain(knproma), &! acc. drainage for HD-model
                               papmecal(knproma)  ! acc. p-e for glaciercalving

    REAL(dp), INTENT(in) :: pros_hd(knproma),   & ! runoff   from *surf* [m]
                            pdrain_hd(knproma), & ! drainage from *surf* [m]
                            palac(knproma),     & ! p - e    from *surf* [m]
                            pdisch(knproma)       ! discharge and calving from
                                                  ! HD and calving model [m/s]
    REAL(dp), INTENT(inout) :: pruntoc(knproma)   ! acc. discharge and calving
                                                  ! for diagnostics kg/(m**2*s)

    ! local scalars

    REAL(dp) ::  zrmean

    ! set accumulated runoff variables zero after HD/coupling time step
    ! (i.e. l_gethd=.true.)
!!!    puthd = io_time_event(1,'days','exact',-delta_time) ! hd model calling frequency
!!!    gethd = io_time_event(1,'days','exact',0)           ! step after puthd

    IF (l_gethd) THEN
      paros(:)    = 0.0_dp
      padrain(:)  = 0.0_dp
      papmecal(:) = 0.0_dp
    END IF

    ! accumulate variables for HD-model

    paros(:)    = paros(:)+pros_hd(:)
    padrain(:)  = padrain(:)+pdrain_hd(:)
    papmecal(:) = papmecal(:)+palac(:)

    ! make means before transfering to HD-Model [m/s]

    IF (l_puthd) THEN
      zrmean = get_interval_seconds(ev_puthd)
      IF (zrmean > 0.0_dp) zrmean = 1.0_dp/zrmean
      paros(:)    = paros(:)*zrmean
      padrain(:)  = padrain(:)*zrmean
      papmecal(:) = papmecal(:)*zrmean
    END IF

    ! accumulate discharge for diagnostics
    ! and convert from m/s to kg/(m**2*s)

    pruntoc(:) = pruntoc(:)+pdisch(:)*rhoh2o*delta_time

  END SUBROUTINE hydrology_collect

  SUBROUTINE hydrology_echam (field_in, field_out, lhd_que)

    !*************************************************************************
    !
    ! **** This program interpolates data from Gaussian grids to a half
    !    degree grid
    !
    !  Programmierung und Entwicklung: Uwe Schulzweida (echamto30min)
    !  Modified to Subroutine by Stefan Hagemann -- September 1995
    !
    ! ***** Version 1.1 -- Dezember 1995
    !          Instead of longitude centred coordinate array trcode, now
    !          the 0.5 degree coordinate field field_out is given back to the calling
    !          routine which has a perfect boundary with the Northpole/dateline
    !          --> Origin has centre coordinate 89.75 N, -179.75 W
    !
    ! ***** Version 2.0 -- November 1999
    !   Since the input array to be transformed is only passed to ECHREAD
    !   but not read in ECHREAD itself, in ECHREAD only
    !   the transformation of the input array to 0.5 degree is done.
    !   Therefor, new calling parameter and routine names are set/given.
    !          alt: SUBROUTINE echread(luinp, ihead, field_out, istep, lhd_que)
    !          neu: SUBROUTINE hydrology_model(field_in, field_out, lhd_que)
    !
    ! ****** List of variables
    !
    !  field_out = Interpolated, transposed Array
    !  lhd_que   = Log-output switch  ( 0 = No Log-Output )

    REAL(dp), INTENT(in) :: field_in(nlon,ngl)
    REAL(dp), INTENT(out) :: field_out(grid_hd%nlon,grid_hd%nlat)
    LOGICAL, INTENT(in) :: lhd_que

    REAL(dp), ALLOCATABLE :: acode(:,:)
    REAL(dp), ALLOCATABLE :: axlon(:)
    REAL(dp), ALLOCATABLE :: axlat(:)
    REAL(dp), ALLOCATABLE :: xr(:), yr(:)

    INTEGER :: jlat, jlon

    ALLOCATE (acode(nlon + 1,ngl + 2))
    ALLOCATE (axlon(nlon + 1))
    ALLOCATE (axlat(ngl + 2))
    ALLOCATE (xr(grid_hd%nlon))
    ALLOCATE (yr(grid_hd%nlat))

    CALL intpol_coord_axis_setup(axlon, axlat, xr, yr)

    IF (lhd_que) THEN
      WRITE(message_text,*) philat(1), philat(2), philat(ngl)
      CALL message('hydrology_echam', message_text)
      WRITE(message_text,*) axlat(1), axlat(2), axlat(ngl + 2)
      CALL message('hydrology_echam', message_text)

      WRITE(message_text,*) philon(1), philon(2), philon(nlon)
      CALL message('hydrology_echam', message_text)
      WRITE(message_text,*) axlon(1), axlon(2), axlon(nlon + 1)
      CALL message('hydrology_echam', message_text)

      WRITE(message_text,*) xr(1), xr(2), xr(grid_hd%nlon)
      CALL message('hydrology_echam', message_text)

      WRITE(message_text,*) yr(1), yr(2), yr(grid_hd%nlat)
      CALL message('hydrology_echam', message_text)
    END IF

    ! generate a copy of field_in with cyclic extension of longitudes and extra
    ! latitudes in the north and south

    DO jlat = 1, ngl
      DO jlon = 1, nlon
        acode(jlon,jlat+1) = field_in(jlon,jlat)
      ENDDO
    ENDDO

    DO jlat = 2, ngl + 1
      acode(nlon + 1,jlat) = acode(1,jlat)
    ENDDO

    DO jlon = 1, nlon + 1
      acode(jlon,1)      = acode(jlon,2)
      acode(jlon,ngl + 2) = acode(jlon,ngl + 2-1)
    ENDDO

    ! interpolation to output grid
    ! (includes transformation on bounded coordinates)
    CALL intpol_with_mapping(nlon + 1, ngl + 2, acode, axlon, axlat, &
         grid_hd%nlon, grid_hd%nlat, field_out, xr, yr, intpol_mapping)

    DEALLOCATE (acode)
    DEALLOCATE (axlon)
    DEALLOCATE (axlat)
    DEALLOCATE (xr)
    DEALLOCATE (yr)

  END SUBROUTINE hydrology_echam

  SUBROUTINE hydrology_corr(nlon, nlat, fatmos, foclsm, aoarea, philat, &
                            fdat, hd_lsm, hd_area, xresi,                   &
                            oclorg, ocscal, florg, fborg, fscal, lhd_que)


    ! **** Correction of the atmospheric grid -> 0.5 degree transformation
    !
    ! ***** Version 1.0 - November 1999
    ! Programmed and developed by Stefan Hagemann, MPI
    ! Remark: Input data FATMOS (Runoff or Drainage) should have the unit m/s.
    !
    ! ***** Version 1.1 - January 2001
    ! Longitude-Index-Correction
    !
    ! ***** Version 2.1 - January 2001
    ! ECHAM5- Version incl. Gaussian latitudes
    !
    !
    ! **** Global Arrays:
    !
    !  fdat = Data array at 0.5 degree
    !  hd_lsm = Land mask array at 0.5 degree
    !  hd_area(jb) = Array of Gridbox areas, Unit = [m^2]
    !
    ! ***** Parameter and arrays of ECHAM grid
    !
    !    nlon = Longitudes of global atmosphere grid
    !    nlat = Latitudes of global atmosphere grid
    !
    !  oclorg = Longitudinal origin of global atmosphere grid
    !  ocscal = resolution/scale = Width of an atmosphere gridbox in degree
    !
    !  fatmos = atmospheric data array
    !  foclsm = land sea mask on atmosphere grid
    !  aoarea = area per gridbox per latitude of atmospheric grid [m^2]
    !   philat = Gaussian latitude of global atmosphere grid (centre coordinates)
    !
    !   xresi = Residuum (Runoff+Drainage), which results from different
    !           land sea masks of the ECHAM and the HD model grid.
    !           In the latest HD model version, a start value is passed to the
    !           HD model that may include further Residual water terms that
    !           should distributed with the discharge to close the water
    !           balance in the coupled atmosphere ocean system.
    !
    !    lhd_que = Log-Output switch ( 0 = No Log-output to STDOUT)
    !
    INTEGER, INTENT(in) :: nlon, nlat
    ! Input fields AO-Grid

    REAL(dp), INTENT(in) :: foclsm(nlon,nlat), aoarea(nlat)
    REAL(dp), INTENT(in) :: fatmos(nlon,nlat), philat(nlat)
    REAL(dp), INTENT(inout) :: fdat(grid_hd%nlon,grid_hd%nlat)
    REAL(dp), INTENT(in) :: hd_lsm(grid_hd%nlon,grid_hd%nlat)
    REAL(dp), INTENT(in) :: hd_area(grid_hd%nlat)
    REAL(dp), INTENT(inout) :: xresi
    REAL(dp), INTENT(in) :: oclorg, ocscal
    REAL(dp), INTENT(in) :: florg, fborg, fscal
    LOGICAL, INTENT(in) :: lhd_que

    REAL(dp) :: x1, x2
    INTEGER :: jb, jl, jlon, jlat

    REAL(dp) :: fb, fl(grid_hd%nlon)
    INTEGER :: jjlon(grid_hd%nlon)
    LOGICAL :: reassigned

    ! set fdat to zero over HD ocean cells
    WHERE (hd_lsm < 0.5)
       fdat = 0._dp
    END WHERE

    ! OA = Ocean-Atmosphere Grid, HD = 0.5 Grad HD-Model Grid
    DO jl = 1, grid_hd%nlon
      fl(jl) = MOD(REAL(jl, dp) * fscal + florg - 0.5_dp * fscal + 180._dp, &
           fullcirc) - 180._dp
    END DO

    DO jl = 1, grid_hd%nlon
      ! Longitude - OCLORG and FL are gridbbox centres
      jjlon(jl) = INT((MOD(fl(jl) - oclorg + ocscal*0.5_dp - fullcirc, &
           -fullcirc) + fullcirc)/ocscal + 1 + 0.00001_dp)
    END DO

    DO jb = 1, grid_hd%nlat-6
      ! First, compute fb and fl for all jl values
      ! grid box centre:
      fb = fborg-REAL(jb,dp)*fscal+0.5_dp*fscal
      ! Corresponding Index in Ocean Grid
      ! Latitude

!!! without Gauss:    XJLAT = (OCBORG-FB) / OCSCAL+1.  , OCBORG = Borderline
!!!                 JLAT = INT(XJLAT+0.00001)

      ! and compute jlat and jlog for all jl values

      jlat = dec_monotonic_closest_midpoint(philat, fb, aub=90._dp, alb=-90._dp)
#if 0
      IF (jlat < 1 .OR. jlat > nlat) THEN
        WRITE(message_text,*) ' error in jlat=', jlat
        CALL message ('hydrology_corr', message_text)
        jlat = nlat
      ENDIF
#endif
      ! laction: Hd_lsm of points, where fdat still has to be computed
      DO jl = 1, grid_hd%nlon
        ! HD Land but OA Water?
        IF (hd_lsm(jl,jb) > 0.5_dp &
             .AND. foclsm(jjlon(jl), jlat) < 0.5_dp) THEN
          jlon = jjlon(jl)
          ! if not --> NW,NE,SW,SE-Directions
          reassigned = reassign_runoff(nlon, nlat, jlon, jlat, &
               foclsm, fatmos, fdat(jl, jb))
          IF (.NOT. reassigned .AND. lhd_que) THEN
            WRITE(message_text,*) 'no land point found for jl=',jl,  &
                 '  jb=',jb, ': fl=',fl(jl), '  fb=',fb
            CALL message('hydrology_corr', message_text)
          END IF
        END IF
      END DO
    ENDDO

    x1 = 0.0_dp
    x2 = 0.0_dp
    DO jl = 1, nlon
      x1 = x1+SUM(fatmos(jl,:)*foclsm(jl,:)*aoarea(:))
    ENDDO
    DO jl = 1, grid_hd%nlon
      x2 = x2+SUM(fdat(jl,:) * hd_lsm(jl,:) * hd_area(:))
    ENDDO
    xresi = xresi+x1-x2

  END SUBROUTINE hydrology_corr
!
! ******************* Interpolation routines ***************************************

  !> find midpoint is,js for every index pair from field
  !> src(src_i_size, src_j_size) to dest(dest_i_size, dest_j_size)
  !> every element cart_idx_2d(is,js) of mapping(id,jd) later defines
  !> that dest(id,jd) will be computed from
  !> src(is,js), src(is - 1, js), src(is, js - 1), src(is - 1, js - 1)
  SUBROUTINE intpol_compute_mapping(mapping, &
       src_i_size, src_j_size, src_x, src_y, &
       dest_i_size, dest_j_size, dest_x, dest_y)

    INTEGER, INTENT(in) :: src_i_size, src_j_size, dest_i_size, dest_j_size

    REAL(dp), INTENT(in)  :: dest_x(dest_i_size),dest_y(dest_j_size)
    REAL(dp), INTENT(in)  :: src_x(src_i_size),src_y(src_j_size)
    TYPE(cart_idx_2d), INTENT(out) :: mapping(dest_i_size,dest_j_size)

    INTEGER :: js, jd, is, id

    mapping = cart_idx_2d(-1, -1)
    DO js = 2, src_j_size
      DO jd = 1, dest_j_size
        ! if (dest_y(jd) < src_y(js-1) .or. dest_y(jd) > src_y(js)) cycle
        IF (dest_y(jd) < MIN(src_y(js - 1), src_y(js)) .OR.   &
            dest_y(jd) > MAX(src_y(js - 1), src_y(js))) CYCLE
        DO is = 2, src_i_size
          DO id = 1, dest_i_size
            IF(dest_x(id) < src_x(is-1) .OR. dest_x(id) > src_x(is)) CYCLE
            mapping(id, jd) = cart_idx_2d(is, js)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE intpol_compute_mapping

  !> interpolate field src to dest
  SUBROUTINE intpol_with_mapping(src_i_size, src_j_size, src, src_x, src_y, &
       dest_i_size, dest_j_size, dest, dest_x, dest_y, mapping)

    INTEGER, INTENT(in) :: src_i_size, src_j_size, dest_i_size, dest_j_size

    REAL(dp), INTENT(in)  :: dest_x(dest_i_size),dest_y(dest_j_size)
    REAL(dp), INTENT(in)  :: src_x(src_i_size),src_y(src_j_size)
    REAL(dp), INTENT(in)  :: src(src_i_size,src_j_size)
    REAL(dp), INTENT(out) :: dest(dest_i_size,dest_j_size)
    TYPE(cart_idx_2d), INTENT(in) :: mapping(dest_i_size, dest_j_size)
    INTEGER :: js, jd, is, id, idt

    dest(:,:) = 0._dp

    DO jd = 1, dest_j_size
      DO id = 1, dest_i_size
        IF (mapping(id,jd)%ilat /= -1) THEN
          is = mapping(id,jd)%ilon
          js = mapping(id,jd)%ilat
          idt = id + dest_i_size/2 - MERGE(dest_i_size, 0, id > dest_i_size/2)
          dest(idt,jd) = src(is - 1, js - 1) &
               * (dest_x(id) - src_x(is)) * (dest_y(jd) - src_y(js)) &
               / ((src_x(is - 1) - src_x(is)) * (src_y(js - 1) - src_y(js))) &
               + src(is, js - 1) &
               * (dest_x(id) - src_x(is - 1)) * (dest_y(jd) - src_y(js)) &
               / ((src_x(is) - src_x(is - 1)) * (src_y(js - 1) - src_y(js))) &
               + src(is - 1, js) &
               * (dest_x(id) - src_x(is)) * (dest_y(jd) - src_y(js - 1)) &
               / ((src_x(is - 1) - src_x(is)) * (src_y(js) - src_y(js - 1))) &
               + src(is, js) &
               * (dest_x(id) - src_x(is - 1)) * (dest_y(jd) - src_y(js - 1)) &
               / ((src_x(is) - src_x(is - 1)) * (src_y(js) - src_y(js - 1)))
        END IF
      END DO
    END DO
  END SUBROUTINE intpol_with_mapping

  FUNCTION reassign_runoff(nlon, nlat, jlon, jlat, foclsm, fatmos, fdat) &
       RESULT(reassigned)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: nlon, nlat, jlon, jlat
    REAL(dp), INTENT(in) :: foclsm(nlon, nlat), fatmos(nlon, nlat)
    REAL(dp), INTENT(inout) :: fdat
    LOGICAL :: reassigned
    !
    ! ndd = If no land point is found as direct neighbour,
    !       it is searched in NWSE direction until the maximum distance of
    !       NDD Boxes is reached.

    INTEGER, PARAMETER :: ndd = 3

    REAL(dp) :: x1, inverted_neighbour_weight
    INTEGER :: idd

    reassigned = .FALSE.
    ! HD Land but OA Water
    ! Considered neighbour gridboxes in OA grid
    ! N,S,W,E-Directions
    x1 = 0.0_dp
    inverted_neighbour_weight = 0.0_dp
    IF (jlon /= 1) THEN
      IF (foclsm(jlon-1,jlat) > 0.5_dp) THEN
        x1 = fatmos(jlon-1,jlat)
        inverted_neighbour_weight = 1.0_dp
      ENDIF
    ELSE
      IF (foclsm(nlon,jlat) > 0.5_dp) THEN
        x1 = fatmos(nlon,jlat)
        inverted_neighbour_weight = 1.0_dp
      ENDIF
    ENDIF
    IF (jlon /= nlon) THEN
      IF (foclsm(jlon+1,jlat) > 0.5_dp) THEN
        x1 = x1+fatmos(jlon+1,jlat)
        inverted_neighbour_weight = inverted_neighbour_weight + 1.0_dp
      ENDIF
    ELSE
      IF (foclsm(1,jlat) > 0.5_dp) THEN
        x1 = x1+fatmos(1,jlat)
        inverted_neighbour_weight = inverted_neighbour_weight + 1.0_dp
      ENDIF
    ENDIF
    IF (jlat /= 1) THEN
      IF (foclsm(jlon,jlat-1) > 0.5_dp) THEN
        x1 = x1+fatmos(jlon,jlat-1)
        inverted_neighbour_weight = inverted_neighbour_weight + 1.0_dp
      ENDIF
    ENDIF
    IF (jlat /= nlat) THEN
      IF (foclsm(jlon,jlat+1) > 0.5_dp) THEN
        x1 = x1+fatmos(jlon,jlat+1)
        inverted_neighbour_weight = inverted_neighbour_weight + 1.0_dp
      ENDIF
    ENDIF
    ! Land point found?
    IF (inverted_neighbour_weight > 0.5_dp) THEN
      fdat = x1/inverted_neighbour_weight
      reassigned = .TRUE.
      RETURN
    END IF

    IF (jlon /= 1) THEN
      IF (jlat /= 1) THEN
        IF (foclsm(jlon-1,jlat-1) > 0.5_dp) THEN
          x1 = x1+fatmos(jlon-1,jlat-1)
          inverted_neighbour_weight = inverted_neighbour_weight + 1.0_dp
        ENDIF
      ENDIF
      IF (jlat /= nlat) THEN
        IF (foclsm(jlon-1,jlat+1) > 0.5_dp) THEN
          x1 = x1+fatmos(jlon-1,jlat+1)
          inverted_neighbour_weight = inverted_neighbour_weight + 1.0_dp
        ENDIF
      ENDIF
    ELSE
      IF (jlat /= 1) THEN
        IF (foclsm(nlon,jlat-1) > 0.5_dp) THEN
          x1 = x1+fatmos(nlon,jlat-1)
          inverted_neighbour_weight = inverted_neighbour_weight + 1.0_dp
        ENDIF
      ENDIF
      IF (jlat /= nlat) THEN
        IF (foclsm(nlon,jlat+1) > 0.5_dp) THEN
          x1 = x1+fatmos(nlon,jlat+1)
          inverted_neighbour_weight = inverted_neighbour_weight + 1.0_dp
        ENDIF
      ENDIF
    ENDIF
    IF (jlon /= nlon) THEN
      IF (jlat /= 1) THEN
        IF (foclsm(jlon+1,jlat-1) > 0.5_dp) THEN
          x1 = x1+fatmos(jlon+1,jlat-1)
          inverted_neighbour_weight = inverted_neighbour_weight + 1.0_dp
        ENDIF
      ENDIF
      IF (jlat /= nlat) THEN
        IF (foclsm(jlon+1,jlat+1) > 0.5_dp) THEN
          x1 = x1+fatmos(jlon+1,jlat+1)
          inverted_neighbour_weight = inverted_neighbour_weight + 1.0_dp
        ENDIF
      ENDIF
    ELSE
      IF (jlat /= 1) THEN
        IF (foclsm(1,jlat-1) > 0.5_dp) THEN
          x1 = x1+fatmos(1,jlat-1)
          inverted_neighbour_weight = inverted_neighbour_weight + 1.0_dp
        ENDIF
      ENDIF
      IF (jlat /= nlat) THEN
        IF (foclsm(1,jlat+1) > 0.5_dp) THEN
          x1 = x1+fatmos(1,jlat+1)
          inverted_neighbour_weight = inverted_neighbour_weight + 1.0_dp
        ENDIF
      ENDIF
    ENDIF
    IF (inverted_neighbour_weight > 0.5_dp) THEN
      ! Second next points in OA grid in N,S,W,E-Directions
      fdat = x1/inverted_neighbour_weight
      reassigned = .TRUE.
      RETURN
    END IF
    extended_surround_loop: DO idd = 2, ndd
      ! HD Land but OA Water
      x1 = 0.0_dp
      inverted_neighbour_weight = 0.0_dp
      IF (jlon-idd >= 1) THEN
        IF (foclsm(jlon-idd,jlat) > 0.5_dp) THEN
          x1 = fatmos(jlon-idd,jlat)
          inverted_neighbour_weight = 1.0_dp
        ENDIF
      ELSE
        IF (foclsm(nlon+jlon-idd,jlat) > 0.5_dp) THEN
          x1 = fatmos(nlon+jlon-idd,jlat)
          inverted_neighbour_weight = 1.0_dp
        ENDIF
      ENDIF
      IF (jlon+idd <= nlon) THEN
        IF (foclsm(jlon+idd,jlat) > 0.5_dp) THEN
          x1 = x1+fatmos(jlon+idd,jlat)
          inverted_neighbour_weight = inverted_neighbour_weight + 1.0_dp
        ENDIF
      ELSE
        IF (foclsm(jlon+idd-nlon,jlat) > 0.5_dp) THEN
          x1 = x1+fatmos(jlon+idd-nlon,jlat)
          inverted_neighbour_weight = inverted_neighbour_weight + 1.0_dp
        ENDIF
      ENDIF
      IF (jlat-idd >= 1) THEN
        IF (foclsm(jlon,jlat-idd) > 0.5_dp) THEN
          x1 = x1+fatmos(jlon,jlat-idd)
          inverted_neighbour_weight = inverted_neighbour_weight + 1.0_dp
        ENDIF
      ENDIF
      IF (jlat+idd <= nlat) THEN
        IF (foclsm(jlon,jlat+idd) > 0.5_dp) THEN
          x1 = x1+fatmos(jlon,jlat+idd)
          inverted_neighbour_weight = inverted_neighbour_weight + 1.0_dp
        ENDIF
      ENDIF
      ! End of Do (IDD) -Loop for Land Point found
      IF (inverted_neighbour_weight > 0.5_dp) THEN
        fdat = x1/inverted_neighbour_weight
        reassigned = .TRUE.
        EXIT extended_surround_loop
      END IF
    END DO extended_surround_loop
    ! end of HD land, OA water
  END FUNCTION reassign_runoff

  SUBROUTINE prepare_intpol_mapping
    REAL(dp) :: axlon(nlon + 1), axlat(ngl + 2), xr(grid_hd%nlon), yr(grid_hd%nlat)
    CALL intpol_coord_axis_setup(axlon, axlat, xr, yr)
    CALL intpol_compute_mapping(intpol_mapping, &
         nlon + 1, ngl + 2, axlon, axlat, grid_hd%nlon, grid_hd%nlat, xr, yr)
  END SUBROUTINE prepare_intpol_mapping

  SUBROUTINE intpol_coord_axis_setup(axlon, axlat, xr, yr)
    REAL(dp), INTENT(inout) :: axlon(nlon + 1), axlat(ngl + 2), xr(grid_hd%nlon), yr(grid_hd%nlat)

    INTEGER :: j

    ! definition of the echam gp data grid
    axlat(1)       =  90.0_dp
    axlat(2:ngl+1) =  philat(1:ngl)
    axlat(ngl + 2) = -90.0_dp

    axlon(1:nlon)  = philon(:)
    axlon(nlon + 1) = fullcirc

    ! definition of hd data grid
    DO j = 1, grid_hd%nlon
      xr(j) = 0.5_dp*grid_hd%resolution + (j-1)/REAL(grid_hd%nlon,dp)*fullcirc
    END DO

    DO j = 1, grid_hd%nlat
      ! perhaps use this formulation? (tj, 20091009)
      ! yr(j) = -(0.5*hd_scal_lat + fullcirc * hd_scal_lat * ( (j-1) &
      !          / REAL(grid_hd%nlat,dp) - 0.5))
      yr(j) = -(0.5_dp*grid_hd%resolution - 0.5_dp*grid_hd%resolution * fullcirc &
           &    + grid_hd%resolution * fullcirc * REAL(j-1, dp) / REAL(grid_hd%nlat,dp))
    ENDDO
  END SUBROUTINE intpol_coord_axis_setup

  !> interpolate field src to dest
  SUBROUTINE intpol (src_i_size, src_j_size, src, src_x, src_y, &
       dest_i_size, dest_j_size, dest, dest_x, dest_y)

    INTEGER, INTENT(in) :: src_i_size, src_j_size, dest_i_size, dest_j_size

    REAL(dp), INTENT(in)  :: dest_x(dest_i_size), dest_y(dest_j_size)
    REAL(dp), INTENT(in)  :: src_x(src_i_size), src_y(src_j_size)
    REAL(dp), INTENT(in)  :: src(src_i_size,src_j_size)
    REAL(dp), INTENT(out) :: dest(dest_i_size,dest_j_size)

    INTEGER :: irun, js, jd, is, id, idt

    irun = 0
    DO js = 2, src_j_size
      DO jd = 1, dest_j_size
        ! if (dest_y(jd) < src_y(js-1) .or. dest_y(jd) > src_y(js)) cycle
        IF (dest_y(jd) < MIN(src_y(js - 1), src_y(js)) .OR.   &
             dest_y(jd) > MAX(src_y(js - 1), src_y(js))) CYCLE
        irun = irun+1
        DO is = 2, src_i_size
          DO id = 1, dest_i_size
            IF(dest_x(id) < src_x(is-1) .OR. dest_x(id) > src_x(is)) CYCLE
            idt = MOD(id - 1 + dest_i_size/2, dest_i_size) + 1
            dest(idt,jd) = src(is - 1, js - 1) &
                 * (dest_x(id) - src_x(is)) * (dest_y(jd) - src_y(js)) &
                 / ((src_x(is - 1) - src_x(is)) * (src_y(js - 1) - src_y(js))) &
                 + src(is, js - 1) &
                 * (dest_x(id) - src_x(is - 1)) * (dest_y(jd) - src_y(js)) &
                 / ((src_x(is) - src_x(is - 1)) * (src_y(js - 1) - src_y(js))) &
                 + src(is - 1, js) &
                 * (dest_x(id) - src_x(is)) * (dest_y(jd) - src_y(js - 1)) &
                 / ((src_x(is - 1) - src_x(is)) * (src_y(js) - src_y(js - 1))) &
                 + src(is, js) &
                 * (dest_x(id) - src_x(is - 1)) * (dest_y(jd) - src_y(js - 1)) &
                 / ((src_x(is) - src_x(is - 1)) * (src_y(js) - src_y(js - 1)))
          ENDDO
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE intpol

END MODULE mo_jsbach_to_hd
