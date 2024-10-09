! mo_interpol.f90 - Utilites for interpolation and discharge conversion 
! 
! Copyright (C) 2021, Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
! SPDX-License-Identifier: Apache-2.0
! See ./LICENSES/ for license information
!
! Authors: Stefan Hagemann
! Contact: <stefan.hagemann@hereon.de>
!_________________________________________

MODULE mo_interpol

  !
  ! Authors:
  !
  ! S. Hagemann, HZG-IfK, May 2018, original source

  use mo_grid,         ONLY: model, distance
!
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: remapnn, generate_point_mapping, mapping, write_mapping_data, &
            island_mapping, define_separation, flow_separation

  !-----------------------------------------------------------------------------

  TYPE mapping
!
!   *** Maximum allowed distances to next target point
    INTEGER :: nmask = 1       ! No. of masks (1 for INEMOU =1,2 or 2 for INEMOU = 3)
                               ! 1 mask -> only secondary mask is used
    DOUBLE PRECISION :: dist_max=100000.      ! 100 km (used for secondary mask if nmask=2)
    DOUBLE PRECISION :: dist_max_prime=25000. ! 25 km for primary mask only if nmask=2
    DOUBLE PRECISION :: deg_max_bound=0.5     ! 0.5 degree if outside target domain
    LOGICAL :: unique = .FALSE.               ! Unique directions F/T = NO/YES
  END TYPE

  TYPE separation_info
!
!   *** Separation coordinates
    INTEGER :: ilonsep    ! Longitude index of gridbox to be separated       
    INTEGER :: ilatsep    ! Latitude index of gridbox to be separated       
    INTEGER :: ilonmouth  ! Longitude index of mouth from where it is separated       
    INTEGER :: ilatmouth  ! Latitude index of  mouth from where it is separated       
  END TYPE
  INTEGER, PARAMETER :: nsepmax = 5         ! Maximum number of separation cases 
  INTEGER, PUBLIC :: nsep = 0                 ! Actual number of separation cases 

  TYPE(separation_info), DIMENSION(nsepmax) :: separation

  !-----------------------------------------------------------------------------
CONTAINS

!   *******************************************************************************
    SUBROUTINE remapnn(nlon, nlat, xlon, xlat, fdata, xmiss, nx, ny, xcoord, ycoord, fout)
!   *******************************************************************************

      INTEGER, INTENT(in)             :: nlon
      INTEGER, INTENT(in)             :: nlat
      DOUBLE PRECISION, DIMENSION(:), INTENT(in)  :: xlon
      DOUBLE PRECISION, DIMENSION(:), INTENT(in)  :: xlat
      REAL, DIMENSION(:,:), INTENT(in)  :: fdata
      REAL, INTENT(in)                :: xmiss
      INTEGER, INTENT(in)             :: nx
      INTEGER, INTENT(in)             :: ny
      REAL, DIMENSION(:,:), INTENT(in)  :: xcoord
      REAL, DIMENSION(:,:), INTENT(in)  :: ycoord
      REAL, DIMENSION(:,:), INTENT(out)  :: fout
!
      DOUBLE PRECISION, PARAMETER          :: zeps = 1.E-10 
      DOUBLE PRECISION, PARAMETER          :: dismax = 100000.  ! 100 km 
      DOUBLE PRECISION, DIMENSION(nx, ny)  :: fdist
      DOUBLE PRECISION                     :: dmin
      REAL :: xmin, xmax, ymin, ymax
      LOGICAL, DIMENSION(nx, ny)  :: lmask
      INTEGER  :: jl, jb
      INTEGER, DIMENSION(2) :: icoord
!
      lmask(:,:) = .true.
      fout(:,:) = 0.
!
      xmin = MINVAL(xcoord)
      xmax = MAXVAL(xcoord)
      ymin = MINVAL(ycoord)
      ymax = MAXVAL(ycoord)
!
      DO jb=1, nlat
      DO jl=1, nlon
      IF (ABS(fdata(jl,jb)).GT.zeps .AND. ABS(fdata(jl,jb)-xmiss).GT.zeps) THEN
      IF (xlon(jl).GE.xmin .AND. xlon(jl).LE.xmax .AND. &
          xlat(jb).GE.ymin .AND. xlat(jb).LE.ymax ) THEN

         CALL distance(xlon(jl), xlat(jb), nx, ny, lmask,    &
              DBLE(xcoord), DBLE(ycoord), fdist)
         dmin = MINVAL(fdist)
         IF (dmin.LE.dismax) THEN
           icoord = MINLOC(fdist)
           fout(icoord(1), icoord(2)) = fdata(jl,jb)
         ENDIF
      ENDIF
      ENDIF
      ENDDO
      ENDDO

    END SUBROUTINE remapnn
!
! *********************************************************************
  SUBROUTINE generate_point_mapping(mask_src, mask_target, map_char, &
             mask_target_prime, xmiss, ique,  &
             mask_src_mapped, ix_target, iy_target, mask_on_ocean)
! *********************************************************************
!
!   Routine mapping (river mouth) points from a source grid (mask_src%value=1) to
!   the nearest river mouth points on a target ocean grid (mask_target%value=1). 
!   It creates arrays on the source grid with the target lon and lat indices that
!   allocate these nearest mouth points to each source mouth point.
!
    TYPE(model), INTENT(in)    :: mask_src     ! Source mask info, e.g. HD mouths 1, otherwise 0
    TYPE(model), INTENT(in)    :: mask_target  ! Source mask info, e.g.
    TYPE(mapping), INTENT(in)  :: map_char     ! Mapping charactistics

    INTEGER, DIMENSION(mask_target%nlon,mask_target%nlat), &
             INTENT(in) :: mask_target_prime  ! Primary mask: Mouths 1, otherwise 0
    DOUBLE PRECISION, INTENT(in)    :: xmiss    ! Missing Value, -9999. is suitable, depends of definition.
    INTEGER, INTENT(in) :: ique     ! Write some debug statements: 0/1 = No/yes
!
    INTEGER, DIMENSION(mask_src%nlon,mask_src%nlat), &
        INTENT(out) :: mask_src_mapped  ! Mask of source mouths mapped to target mouth
    INTEGER, DIMENSION(mask_src%nlon,mask_src%nlat), &
        INTENT(out) :: ix_target        ! x-Indices of nearest target mouth points
    INTEGER, DIMENSION(mask_src%nlon,mask_src%nlat), &
        INTENT(out) :: iy_target        ! y-Indices of nearest target mouth points
    INTEGER, DIMENSION(mask_target%nlon,mask_target%nlat), &
        INTENT(out) :: mask_on_ocean    ! Mask on ocean target receiving mapping
!
    DOUBLE PRECISION :: zlon, zlat
    LOGICAL, DIMENSION(:,:), ALLOCATABLE :: lmask1  ! Logical array of primary mask
    LOGICAL, DIMENSION(:,:), ALLOCATABLE :: lmask2  ! Logical array of secondary mask
    LOGICAl :: lfound
    INTEGER :: jl, jb, nmes, nmou
    INTEGER :: nomo                       ! Number of HD points with no ocean model mouth found
    DOUBLE PRECISION :: xmax, xmin, ymax, ymin
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: fdist
    DOUBLE PRECISION :: dmin
    INTEGER, DIMENSION(2) :: icoord
!   *** Variables for achieving unique pointers
    INTEGER :: nmix  ! nmix = Maximum number of HD mouths pointing to a single Oceanmodel mouth
    INTEGER jx, jy, ibox, I, nact, il,ib
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: dmix
    INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: icmix
    LOGICAL, DIMENSION(:,:), ALLOCATABLE :: lfree
    DOUBLE PRECISION, PARAMETER :: zeps = 1.E-4
!
    ALLOCATE(fdist(mask_target%nlon,mask_target%nlat))
!
!   *** Domain Boundaries
    xmin = MINVAL(mask_target%xlon)
    xmax = MAXVAL(mask_target%xlon)
    ymin = MINVAL(mask_target%xlat)
    ymax = MAXVAL(mask_target%xlat)
    WRITE(*,*) "Ocean Model-Region: Lon ", xmin,' - ', xmax 
    WRITE(*,*) "                    Lat ", ymin,' - ', ymax 
    WRITE(*,*)        "First - Last Lat ", mask_target%xlat(1,1),' - ', &
                mask_target%xlat(mask_target%nlon,mask_target%nlat) 
!
!   *** Note that (opposite to NEMO) mouth masks defines mouth points as 1.
    ALLOCATE(lmask2(mask_target%nlon,mask_target%nlat))
    lmask2(:,:) = .FALSE.
    WHERE(mask_target%value.EQ.1)        
      lmask2=.true.
    END WHERE
    IF (map_char%nmask.EQ.2) THEN
      ALLOCATE(lmask1(mask_target%nlon,mask_target%nlat))
      lmask1(:,:) = .FALSE.
      WHERE(mask_target_prime.EQ.1)        
        lmask1=.true.
      END WHERE
    ENDIF
!
!   ******** Searching for nearest NEMO mouth point
    ix_target(:,:) = 0
    iy_target(:,:) = 0
    mask_src_mapped(:, :) = NINT(xmiss)
    mask_on_ocean(:,:)=0
    nmes=0
    nmou=0
    nomo = 0
    DO jb = 1, mask_src%nlat
    DO jl = 1, mask_src%nlon
    IF ( ABS(mask_src%value(jl,jb)-1).LE.0.1 ) THEN   ! Mouth point = 1 
      IF (mask_src%xlon(jl,jb).GE.xmin-map_char%deg_max_bound .AND. &
          mask_src%xlon(jl,jb).LE.xmax+map_char%deg_max_bound .AND. &
          mask_src%xlat(jl,jb).GE.ymin-map_char%deg_max_bound .AND. &
          mask_src%xlat(jl,jb).LE.ymax+map_char%deg_max_bound) THEN
        nmou = nmou+1
        mask_src_mapped(jl, jb) = 0
        lfound = .FALSE.
        IF (map_char%nmask.EQ.2) THEN
  	  CALL DISTANCE(mask_src%xlon(jl,jb), mask_src%xlat(jl,jb), &
                        mask_target%nlon, mask_target%nlat, lmask1,  &
                        mask_target%xlon, mask_target%xlat, fdist)
          icoord = MINLOC(fdist, MASK=lmask1)
          dmin = MINVAL(fdist, MASK=lmask1)
          IF (dmin.LE.map_char%dist_max_prime) lfound = .TRUE.
          IF (icoord(2).GT.730 .AND. icoord(1).GT.700) THEN
            WRITE(*, '(A,2(X,I4),A, 2(X,F7.2), A,X,F6.2,A,2(X,I4) )') &
                'M1: - Lon/Lat: ', jl, jb, ' = ', mask_src%xlon(jl,jb), mask_src%xlat(jl,jb), &
                ' Distance HD-Oceanmodel mouth = ', dmin/1000., ' km at ', icoord
            WRITE(*,*) '     -> Lon = ', mask_target%xlon(icoord(1), icoord(2)), &
                 ' Lat = ', mask_target%xlat(icoord(1),icoord(2))
          ENDIF
        ENDIF
        IF (.NOT. lfound) THEN
          CALL DISTANCE(mask_src%xlon(jl,jb), mask_src%xlat(jl,jb), &
                        mask_target%nlon, mask_target%nlat, lmask2,  &
                        mask_target%xlon, mask_target%xlat, fdist)
          icoord = MINLOC(fdist, MASK=lmask2)
          dmin = MINVAL(fdist, MASK=lmask2)
          IF (dmin.LE.map_char%dist_max) lfound = .TRUE.
          IF (ique.ne.0 .AND. icoord(2).GT.730 .AND. icoord(1).GT.700) THEN
            WRITE(*, '(A,2(X,I4),A, 2(X,F7.2), A,X,F6.2,A,2(X,I4) )') &
                'M2 - Lon/Lat: ', jl, jb, ' = ', mask_src%xlon(jl,jb), mask_src%xlat(jl,jb),  &
                ' Distance HD-Oceanmodel mouth = ', dmin/1000., ' km at ', icoord
            WRITE(*,*) '     -> Lon = ', mask_target%xlon(icoord(1), icoord(2)), &
                ' Lat = ', mask_target%xlat(icoord(1),icoord(2))
          ENDIF
        ENDIF
!
        IF (lfound) THEN
          mask_src_mapped(jl, jb) = 1
          nmes = nmes+1
          ix_target(jl,jb) = icoord(1)
          iy_target(jl,jb) = icoord(2)
          mask_on_ocean(icoord(1), icoord(2)) = mask_on_ocean(icoord(1), icoord(2)) + 1
          IF (ique.EQ.1) THEN
            WRITE(*,*) nmes, '. Distance Source mouth -> Target mouth = ', &
                       dmin/1000., ' km at ',  icoord
            WRITE(*,*) '     -> Lon = ', mask_target%xlon(icoord(1), icoord(2)), &
                       ' Lat = ', mask_target%xlat(icoord(1),icoord(2))
          ENDIF
        ELSE
          nomo = nomo + 1
        ENDIF
      ENDIF
    ENDIF
    ENDDO
    ENDDO
    WRITE(*,*) nmou, ' Source (e.g. HD) Points in ',TRIM(mask_target%name), &
        ' region --> Converted to ', nmes, ' ', TRIM(mask_target%name),' mouth points'
    nmix = MAXVAL(mask_on_ocean)
    WRITE(*,*) 'Maximum number of source mouths pointing towards a single ', &
               TRIM(mask_target%name),' mouth: ', nmix
    WRITE(*,*) 'No ', TRIM(mask_target%name), ' mouth was found for ', nomo, ' source (e.g. HD) mouths'
!
!   *** Indicate source mouths that point towards a Ocean model mouth with more than one input  
    DO jb = 1, mask_src%nlat
    DO jl = 1, mask_src%nlon
    IF ( mask_src_mapped(jl,jb).EQ.1 ) THEN
      IF (mask_on_ocean(ix_target(jl,jb), iy_target(jl,jb)).GT.1)   &
          mask_src_mapped(jl, jb) = mask_on_ocean(ix_target(jl,jb), iy_target(jl,jb))
    ENDIF
    ENDDO
    ENDDO
!
    IF (map_char%unique) THEN
      IF (nmix.GT.2) THEN
        WRITE(*,*) "****** ERROR - Cases with more more than 2 HD mouths: up to ", nmix
        WRITE(*,*) "       pointing to the ",TRIM(mask_target%name)," NEMo box are NOT implemented"
        STOP       " --> Program stops right here!!!!"
      ENDIF
      ALLOCATE( dmix(nmix, nmix) )
      ALLOCATE( icmix(nmix, nmix, 2) )
      ALLOCATE( lfree(mask_target%nlon,mask_target%nlat) )
      WHERE (mask_on_ocean(:,:).EQ.1) lmask2(:,:) = .FALSE.
      DO jb = 1, mask_src%nlat
      DO jl = 1, mask_src%nlon
        IF ( mask_src_mapped(jl,jb).EQ.2 ) THEN
          jx = ix_target(jl,jb)
          jy = iy_target(jl,jb)
          nact = mask_on_ocean(jx,jy)
          zlon = mask_src%xlon(jl,jb)
          zlat = mask_src%xlat(jl,jb)
          lfree(:,:) = lmask2(:,:)
          I=0
          WRITE(*, '(A,2(X,I4),A,2(X,F7.2) ,A,2(X,F7.2),A, 2(X,I4) )') &
                   'Mixpoint - Lon/Lat: ', jl,jb, ' = ', zlon, zlat, ' pointing towards ', &
                    mask_target%xlon(jx,jy), mask_target%xlat(jx,jy), ' --> jx, jy: ', jx, jy
!
!         *** Distances 1, nact
          DO WHILE(I.LT.nact)
            I=I+1
	    CALL DISTANCE(zlon, zlat, mask_target%nlon, mask_target%nlat, lfree,  &
               mask_target%xlon, mask_target%xlat, fdist)
            icmix(1,I,:) = MINLOC(fdist, MASK=lfree)
            dmix(1,I) = MINVAL(fdist, MASK=lfree)
            lfree(icmix(1,I,1),icmix(1,I,2)) = .FALSE.
          ENDDO
         
          IF (ique.EQ.1) WRITE(*,*) 'Dist1 ', dmix(1,1:nact) 
!
!         *** HD mouths poiting to the same Ocean model mouth
          ibox=2
          DO ib = 1, mask_src%nlat
          DO il = 1, mask_src%nlon
          IF (mask_src_mapped(il,ib).EQ.2 .AND. (il.NE.jl .OR. ib.NE.jb) .AND.  &
             ix_target(il,ib).EQ.jx .AND. iy_target(il,ib).EQ.jy ) THEN

            IF (ique.EQ.1) WRITE(*,*) 'il, ib ', il, ib,  ' -> ', &
                               ix_target(il,ib), iy_target(il,ib)
            zlon = mask_src%xlon(il,ib)
            zlat = mask_src%xlat(il,ib)
            lfree(:,:) = lmask2(:,:)

            I=0
!
!           *** Distances ibox, nact
            DO WHILE(I.LT.nact)
              I=I+1
	      CALL DISTANCE(zlon, zlat, mask_target%nlon, mask_target%nlat, lfree,  &
                 mask_target%xlon, mask_target%xlat, fdist)
              icmix(ibox,I,:) = MINLOC(fdist, MASK=lfree)
              dmix(ibox,I) = MINVAL(fdist, MASK=lfree)
              lfree(icmix(ibox,I,1),icmix(ibox,I,2)) = .FALSE.
            ENDDO

            IF (ique.EQ.1) WRITE(*,*) 'Dist2 ', dmix(2,1:nact) 
!
!           *** Only nmix = nact = 2 can be dealt with
            IF (dmix(1,1).LE.dmix(2,1)) THEN           ! Box 1 closer than 2
              mask_src_mapped(jl,jb) = 1           
              IF (dmix(2,2).LE.map_char%dist_max) THEN
                ix_target(il,ib) = icmix(2,2,1)
                iy_target(il,ib) = icmix(2,2,2)
                mask_on_ocean(icmix(2,2,1), icmix(2,2,2)) = 1
                lmask2(icmix(2,2,1), icmix(2,2,2)) = .FALSE.
!!                IF (ique.EQ.1) THEN
                WRITE(*,*) 'New Distance Source-Target mouth = ', &
                           dmix(2,2)/1000., ' km at ',  icmix(2,2,:)
                IF (ique.EQ.1) WRITE(*,*) '     -> Lon = ', mask_target%xlon(icmix(2,2,1), &
                       icmix(2,2,2)),  ' Lat = ', mask_target%xlat(icmix(2,2,1), icmix(2,2,2))
!!                ENDIF
                mask_src_mapped(il,ib) = 1
              ELSE
                mask_src_mapped(il,ib) = 0
                ix_target(il,ib) = 0
                iy_target(il,ib) = 0
                WRITE(*,*) 'Second closest target mouth is out of range: ', &
                           dmix(2,2)/1000., ' km at ',  icmix(2,2,:)
              ENDIF
            ELSE                                       ! Box 2 closer than 1
              mask_src_mapped(il,ib) = 1           
              IF (dmix(1,2).LE.map_char%dist_max) THEN
                ix_target(jl,jb) = icmix(1,2,1)
                iy_target(jl,jb) = icmix(1,2,2)
                mask_on_ocean(icmix(1,2,1), icmix(1,2,2)) = 1
                lmask2(icmix(1,2,1), icmix(1,2,2)) = .FALSE.
!!                IF (ique.EQ.1) THEN
                WRITE(*,*) 'New Distance Source-targetNEMO mouth = ', &
                           dmix(1,2)/1000., ' km at ',  icmix(1,2,:)
                IF (ique.EQ.1) WRITE(*,*) '     -> Lon = ', mask_target%xlon(icmix(1,2,1), &
                     icmix(1,2,2)),  ' Lat = ', mask_target%xlat(icmix(1,2,1), icmix(1,2,2))
!!                ENDIF
                mask_src_mapped(jl,jb) = 1
              ELSE
                mask_src_mapped(jl,jb) = 0
                ix_target(jl,jb) = 0
                iy_target(jl,jb) = 0
                WRITE(*,*) 'Second closest target mouth is out of range: ', &
                           dmix(1,2)/1000., ' km at ',  icmix(1,2,:)
              ENDIF
            ENDIF                           
!
!           *** Mixed problem solved for jx,jy
            mask_on_ocean(jx, jy) = 1
            lmask2(jx, jy) = .FALSE.
            EXIT
          ENDIF
          ENDDO
          ENDDO
        ENDIF
      ENDDO
      ENDDO
      DEALLOCATE(dmix)
      DEALLOCATE(icmix)
      DEALLOCATE(lfree)
    ENDIF

    DEALLOCATE(fdist)
    IF (map_char%nmask.EQ.2) DEALLOCATE(lmask1)
    DEALLOCATE(lmask2)

  END SUBROUTINE generate_point_mapping

! *********************************************************************
  SUBROUTINE island_mapping(mask_src, mask_target, map_char, ique,  &
              mask_src_mapped, ix_target, iy_target, mask_on_ocean)
! *********************************************************************
!
!   Routine that maps river mouths from smaller islands to ocean points, which were not 
!   mapped before. 
!   If the ocean model resolution is coarser than the source resolution, some islands on the source 
!   grid may not have coastal ocean points nearby. 
!
    TYPE(model), INTENT(in)    :: mask_src     ! Source mask info, e.g. HD mouths 1, otherwise 0
    TYPE(model), INTENT(in)    :: mask_target  ! Source mask info, e.g.
    TYPE(mapping), INTENT(in)  :: map_char     ! Mapping charactistics

    INTEGER, INTENT(in) :: ique     ! Write some debug statements: 0/1 = No/yes
!
    INTEGER, DIMENSION(mask_src%nlon,mask_src%nlat), &
        INTENT(inout) :: mask_src_mapped  ! Mask of source mouths mapped to target mouth
    INTEGER, DIMENSION(mask_src%nlon,mask_src%nlat), &
        INTENT(inout) :: ix_target        ! x-Indices of nearest target mouth points
    INTEGER, DIMENSION(mask_src%nlon,mask_src%nlat), &
        INTENT(inout) :: iy_target        ! y-Indices of nearest target mouth points
    INTEGER, DIMENSION(mask_target%nlon,mask_target%nlat), &
        INTENT(inout) :: mask_on_ocean    ! Mask on ocean target receiving mapping
!
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: fdist
    DOUBLE PRECISION :: dmin
    INTEGER, DIMENSION(2) :: icoord
!
    INTEGER :: jl, jb, jx, jy, nisland
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: target_island
    LOGICAL, DIMENSION(:,:), ALLOCATABLE :: lmask  ! Logical array
    DOUBLE PRECISION, PARAMETER :: zeps = 1.E-4
!
    ALLOCATE(target_island(mask_target%nlon,mask_target%nlat))
    target_island(:,:) = 0
    ALLOCATE(fdist(mask_target%nlon,mask_target%nlat))
    ALLOCATE( lmask(mask_target%nlon,mask_target%nlat) )
    WHERE(mask_target%value.EQ.0)        
      lmask=.true.
    END WHERE

    nisland = 0
    DO jb = 1, mask_src%nlat
    DO jl = 1, mask_src%nlon
      IF ( mask_src_mapped(jl,jb).EQ.0 ) THEN
        CALL DISTANCE(mask_src%xlon(jl,jb), mask_src%xlat(jl,jb), &
                        mask_target%nlon, mask_target%nlat, lmask,  &
                        mask_target%xlon, mask_target%xlat, fdist)
        icoord = MINLOC(fdist, MASK=lmask)
        dmin = MINVAL(fdist, MASK=lmask)
        IF (dmin.LE.map_char%dist_max) THEN
          ix_target(jl,jb) = icoord(1)
          iy_target(jl,jb) = icoord(2)
          mask_on_ocean(icoord(1), icoord(2)) = mask_on_ocean(icoord(1), icoord(2)) + 1
          target_island(icoord(1), icoord(2)) = 1
          mask_src_mapped(jl, jb) = mask_on_ocean(icoord(1), icoord(2))
          nisland = nisland + 1
          IF (ique.GT.1) THEN
            WRITE(*, '(I3, A,2(X,I4),A, 2(X,F7.2), A,X,F6.2,A,2(X,I4) )') &
                   nisland, '. Island - Lon/Lat: ', jl, jb, ' = ', mask_src%xlon(jl,jb), mask_src%xlat(jl,jb),  &
                   ' Distance HD-Oceanmodel mouth = ', dmin/1000., ' km at ', icoord
            WRITE(*,*) '     -> Lon = ', mask_target%xlon(icoord(1), icoord(2)), &
                   ' Lat = ', mask_target%xlat(icoord(1),icoord(2))
          ENDIF
        ENDIF
      ENDIF
    ENDDO
    ENDDO
    WRITE(*,*) nisland, ' source island points in ',TRIM(mask_target%name), &
        ' region --> Converted to ', SUM(target_island), ' ', TRIM(mask_target%name),' ocean points'
   
    DEALLOCATE(fdist)
    DEALLOCATE(lmask)

  END SUBROUTINE island_mapping

! *********************************************************************
  SUBROUTINE write_mapping_data(dnout, mask_src, mask_target, imode, &
             mask_src_mapped, ix_target, iy_target, mask_on_ocean)
! *********************************************************************
!
!   Routine that writes the mapping arrays created by routine generate_point_mapping) 
!   into dnout. These arrays map the source grid mouths to the target lon and lat 
!   indices the nearest target mouth points.
!
    use netcdf
    use mo_grid,         ONLY: model
	  
    CHARACTER (LEN=*), INTENT(in)  :: dnout       ! Output file name
    TYPE(model), INTENT(in)        :: mask_src    ! Source info & coordinates
    TYPE(model), INTENT(in)        :: mask_target ! Target info, e.g.
    INTEGER :: imode             !  --> Method of using potential mouth masks
                   ! 1    Use existing mask with potential mouth points on ocean grid
                   ! 2    Generate mask of coastal ocean points from ocean sea mask (Def.)
                   ! 3    Combine methods 1 and 2
                   ! 4    as 3, but with both masks prescribed
!
    INTEGER, DIMENSION(mask_src%nlon,mask_src%nlat), &
        INTENT(in) :: mask_src_mapped  ! Mask of source mouths mapped to target mouth
    INTEGER, DIMENSION(mask_src%nlon,mask_src%nlat), &
        INTENT(in) :: ix_target        ! x-Indices of nearest target mouth points
    INTEGER, DIMENSION(mask_src%nlon,mask_src%nlat), &
        INTENT(in) :: iy_target        ! y-Indices of nearest target mouth points
    INTEGER, DIMENSION(mask_target%nlon,mask_target%nlat), &
        INTENT(in) :: mask_on_ocean    ! Mask on ocean target receiving mapping
!
!   *** NETCDF variables
    INTEGER :: ierr, ncid, varid, dimids_src(2), dimids_target(2)
    CHARACTER (LEN=120) :: clong 

    CHARACTER(8) :: cdate     ! Date in reality = date of creating output file

    ierr = nf90_create(dnout, NF90_CLOBBER, ncid)
    ierr = nf90_def_dim(ncid, 'lon', mask_src%nlon, dimids_src(1))
    ierr = nf90_def_dim(ncid, 'lat', mask_src%nlat, dimids_src(2))

    ierr = nf90_def_var(ncid, 'FMOU_HD_TO_NEMO', NF90_INT, dimids_src, varid)
    ierr = nf90_put_att(ncid,varid,'code', 731)
    ierr = nf90_put_att(ncid,varid,'units','[]')
    clong = TRIM(mask_src%name) // ' river mouths pointing towards ' // TRIM(mask_target%name) // ' ocean grid'
    ierr = nf90_put_att(ncid,varid,'long_name',TRIM(clong))
    ierr = nf90_put_att(ncid,varid,'standard_name','FMOU_HD_TO_NEMO')
    ierr = nf90_put_att(ncid,varid,'coordinates','lon lat')
    ierr = nf90_put_att(ncid,varid,'missing_value',-9999.)
    ierr = nf90_enddef(ncid)
    ierr = nf90_put_var(ncid, varid, mask_src_mapped)
!
    ierr = nf90_redef(ncid)
    ierr = nf90_def_var(ncid, 'INDEXX', NF90_INT, dimids_src, varid)
    ierr = nf90_put_att(ncid,varid,'code', 732)
    ierr = nf90_put_att(ncid,varid,'units','[]')
    clong = TRIM(mask_target%name) // ' ocean model longitude index of ' // TRIM(mask_src%name) // ' river mouth'
    ierr = nf90_put_att(ncid,varid,'long_name', TRIM(clong))
    ierr = nf90_put_att(ncid,varid,'standard_name','INDEXX')
    ierr = nf90_put_att(ncid,varid,'coordinates','lon lat')
    ierr = nf90_enddef(ncid)
    ierr = nf90_put_var(ncid, varid, ix_target)

    ierr = nf90_redef(ncid)
    ierr = nf90_def_var(ncid, 'INDEXY', NF90_INT, dimids_src, varid)
    ierr = nf90_put_att(ncid,varid,'code', 733)
    ierr = nf90_put_att(ncid,varid,'units','[]')
    clong = TRIM(mask_target%name) // ' ocean model latitude index of ' // TRIM(mask_src%name) // ' river mouth'
    ierr = nf90_put_att(ncid,varid,'long_name', TRIM(clong))
    ierr = nf90_put_att(ncid,varid,'standard_name','INDEXY')
    ierr = nf90_put_att(ncid,varid,'coordinates','lon lat')
    ierr = nf90_enddef(ncid)
    ierr = nf90_put_var(ncid, varid, iy_target)

    IF (mask_src%idim.EQ.1) THEN
      ierr = nf90_redef(ncid)
      ierr = nf90_def_var(ncid, 'lon', NF90_DOUBLE, dimids_src(1), varid)
      ierr = nf90_put_att(ncid,varid,'units','degrees_east')
      ierr = nf90_put_att(ncid,varid,'long_name','longitude')
      ierr = nf90_put_att(ncid,varid,'standard_name','longitude')
      ierr = nf90_enddef(ncid)
      ierr = nf90_put_var(ncid, varid, mask_src%xlon(:,1))

      ierr = nf90_redef(ncid)
      ierr = nf90_def_var(ncid, 'lat', NF90_DOUBLE, dimids_src(2), varid)
      ierr = nf90_put_att(ncid,varid,'units','degrees_north')
      ierr = nf90_put_att(ncid,varid,'long_name','latitude')
      ierr = nf90_put_att(ncid,varid,'standard_name','latitude')
      ierr = nf90_enddef(ncid)
      ierr = nf90_put_var(ncid, varid, mask_src%xlat(1,:))     !TODO: will not work for ICON
    ELSE
      ierr = nf90_redef(ncid)
      ierr = nf90_def_var(ncid, 'lon', NF90_DOUBLE, dimids_src, varid)
      ierr = nf90_put_att(ncid,varid,'units','degrees_east')
      ierr = nf90_put_att(ncid,varid,'long_name','longitude')
      ierr = nf90_put_att(ncid,varid,'standard_name','longitude')
      ierr = nf90_enddef(ncid)
      ierr = nf90_put_var(ncid, varid, mask_src%xlon)

      ierr = nf90_redef(ncid)
      ierr = nf90_def_var(ncid, 'lat', NF90_DOUBLE, dimids_src, varid)
      ierr = nf90_put_att(ncid,varid,'units','degrees_north')
      ierr = nf90_put_att(ncid,varid,'long_name','latitude')
      ierr = nf90_put_att(ncid,varid,'standard_name','latitude')
      ierr = nf90_enddef(ncid)
      ierr = nf90_put_var(ncid, varid, mask_src%xlat)
    ENDIF
    WRITE(*,*) "River-Mouths pointing towards ocean grid were written: ", TRIM(dnout)

    ierr = nf90_redef(ncid)
    ierr = nf90_def_dim(ncid, 'x', mask_target%nlon, dimids_target(1))
    ierr = nf90_def_dim(ncid, 'y', mask_target%nlat, dimids_target(2))

    ierr = nf90_def_var(ncid, 'FMOU_HD_ON_NEMO', NF90_INT, dimids_target, varid)
    ierr = nf90_put_att(ncid,varid,'code', 741)
    ierr = nf90_put_att(ncid,varid,'units','[]')
    clong = TRIM(mask_src%name) // ' river mouth points on ' // TRIM(mask_target%name) // ' ocean model grid'
    ierr = nf90_put_att(ncid,varid,'long_name', TRIM(clong))
    ierr = nf90_put_att(ncid,varid,'standard_name','FMOU_HD_ON_NEMO')
    ierr = nf90_put_att(ncid,varid,'coordinates','olon olat')
    ierr = nf90_enddef(ncid)
    ierr = nf90_put_var(ncid, varid, mask_on_ocean)

    ierr = nf90_redef(ncid)
    ierr = nf90_def_var(ncid, 'olon', NF90_DOUBLE, dimids_target, varid)
    ierr = nf90_put_att(ncid,varid,'units','degrees_E')
    ierr = nf90_put_att(ncid,varid,'long_name','longitude')
    ierr = nf90_put_att(ncid,varid,'standard_name','longitude')
    ierr = nf90_enddef(ncid)
    ierr = nf90_put_var(ncid, varid, mask_target%xlon)

    ierr = nf90_redef(ncid)
    ierr = nf90_def_var(ncid, 'olat', NF90_DOUBLE, dimids_target, varid)
    ierr = nf90_put_att(ncid,varid,'units','degrees_N')
    ierr = nf90_put_att(ncid,varid,'long_name','latitude')
    ierr = nf90_put_att(ncid,varid,'standard_name','latitude')
    ierr = nf90_enddef(ncid)
    ierr = nf90_put_var(ncid, varid, mask_target%xlat)
!
!   *** Write global attributes
    ierr = nf90_redef(ncid)
    clong = "Mapping data from " // TRIM(mask_src%name) // " mouth points to " // &
             TRIM(mask_target%name) // " coastal points"
    ierr = nf90_put_att(ncid, NF90_GLOBAL, 'Version', clong)
    IF (imode.EQ.1) THEN
      clong = "Potential " // TRIM(mask_target%name) // " mouth (=coastal) points derived from " // &
               "existing mask, e.g. EHYPE."
    ELSE IF (imode.EQ.2) THEN
      clong = "Potential " // TRIM(mask_target%name) // " mouth (=coastal) points derived from " // &
               TRIM(mask_target%name) // " sea mask."
    ELSE IF (imode.EQ.3) THEN
      clong = "Potential " // TRIM(mask_target%name) // " mouth (=coastal) points derived from " // &
               TRIM(mask_target%name) // " sea mask and existing prescribed mask."
    ELSE IF (imode.EQ.4) THEN
      clong = "Potential " // TRIM(mask_target%name) // " mouth (=coastal) points derived from " // &
               TRIM(mask_target%name) // " prescribed primary and secondary masks."
    ENDIF
    CALL DATE_AND_TIME(cdate)          ! returns YYYYMMDD

    ierr = nf90_put_att(ncid, NF90_GLOBAL, 'Conventions',  'CF-1.6')
    ierr = nf90_put_att(ncid, NF90_GLOBAL, 'source',       'HD model Vs. 5, doi:10.5281/zenodo.4893099')
    ierr = nf90_put_att(ncid, NF90_GLOBAL, 'title',        'HD model coupling file')
    ierr = nf90_put_att(ncid, NF90_GLOBAL, 'licence',      'CC-BY 4.0')
    ierr = nf90_put_att(ncid, NF90_GLOBAL, 'creation_date', cdate(7:8)//'.'//cdate(5:6)//'.'//cdate(1:4) )
    ierr = nf90_put_att(ncid, NF90_GLOBAL, 'originator',   'Stefan Hagemann, Hereon' )
    ierr = nf90_put_att(ncid, NF90_GLOBAL, 'contact',      'stefan.hagemann@hereon.de' )

    ierr = nf90_put_att(ncid, NF90_GLOBAL, 'Comment1', TRIM(clong))
    ierr = nf90_put_att(ncid, NF90_GLOBAL, 'Comment2', "FMOU_HD_TO_NEMO may comprise values of")
    ierr = nf90_put_att(ncid, NF90_GLOBAL, 'Comment3', &
           "  0  Mouth in ocean domain but no coastal mouth point nearby.")
    ierr = nf90_put_att(ncid, NF90_GLOBAL, 'Comment4', &
           "  1  Mouth pointing to coastal ocean mouth point with no other inflow source.")
    ierr = nf90_put_att(ncid, NF90_GLOBAL, 'Comment5', &
           " >2  Mouth pointing to coastal ocean mouth point with other inflow sources.")

    ierr = nf90_enddef(ncid)

    ierr = nf90_close(ncid)
    WRITE(*,*) "River-Mouths on ocean model grid were written: ", TRIM(dnout)

  END SUBROUTINE write_mapping_data

! *********************************************************************
  SUBROUTINE define_separation(ilonsep, ilatsep, ilonmouth, ilatmouth)
! *********************************************************************
!
! *** Define flow separation: take water from defined grid boxes (the flows to be separated and 
! *** considered as new mouth boxes in convert_discharge) and subtract the
! *** discharge amounts from the corresponding mouth boxes of the respective parent catchment. 
    INTEGER, INTENT(in) :: ilonsep, ilatsep, ilonmouth, ilatmouth

    nsep = nsep + 1
    separation(nsep)%ilonsep = ilonsep
    separation(nsep)%ilatsep = ilatsep
    separation(nsep)%ilonmouth = ilonmouth
    separation(nsep)%ilatmouth = ilatmouth
!
  END SUBROUTINE define_separation

! *********************************************************************
  SUBROUTINE flow_separation(NL, NB, friv_src)
! *********************************************************************
!
! *** Flow separation: take water from defined grid boxes (the flows to be separated and 
! *** considered as new mouth boxes in convert_discharge) and subtract the
! *** discharge amounts from the corresponding mouth boxes of the respective parent catchment. 
!
    INTEGER, INTENT(in) :: NL, NB
    REAL, DIMENSION(NL,NB), INTENT(inout) :: friv_src
    REAL, DIMENSION(NL,NB) :: fdum
    INTEGER :: i

    fdum = friv_src
    DO i=1, nsep
      friv_src(separation(i)%ilonmouth, separation(i)%ilatmouth) = &
         friv_src(separation(i)%ilonmouth, separation(i)%ilatmouth) - fdum(separation(i)%ilonsep, separation(i)%ilatsep)
    ENDDO

  END SUBROUTINE flow_separation
!
END MODULE mo_interpol
