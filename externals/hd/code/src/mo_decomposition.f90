! mo_decomposition.f90 - Definition of data types for data decomposition in parallel mode
!
! Copyright (C) 2014, MPI-M
! SPDX-License-Identifier: BSD-3-Clause
! See ./LICENSES/ for license information
!_________________________________________

MODULE mo_decomposition
  !
  ! This module defines the data types holding meta information
  !   for data decomposition in parallel mode
  !
  ! Authors:
  !
  ! A. Rhodin, MPI, August 1999, original source
  ! A. Rhodin, DWD/MPI, August 2001, ffsl decomposition 
  ! T. Diehl,  DKRZ, September 2001, spitfire decomposition
  ! L. Kornblueh, merging of ffsl, spitfire 
  ! A. Rhodin, DWD/MPI, May 2002, blocking (nproma)
  !
  !  This routine originates (year 2014) from MPI-ESM, the Earth System Model of the 
  !  Max Planck Institute for Meteorology (Mauritsen et al. 2019). 
  !  Reference: Mauritsen, T., et al. (2019) Developments in the MPI-M Earth System Model 
  !  version 1.2 (MPI-ESM1.2) and its response to increasing CO2. J. Adv. Model. Earth Syst., 11, 
  !  doi: 10.1029/2018MS001400.

  USE mo_exception,   ONLY: finish, message, message_text
  USE mo_util_string, ONLY: separator
  USE mo_mpi,         ONLY: p_nprocs, p_set_communicator, p_pe, p_io

  IMPLICIT NONE

  PRIVATE
  !
  ! Data type declaration
  !
  PUBLIC :: pe_decomposed        ! data type with decomposition of a single PE 
  PUBLIC :: default_decomposed
  PUBLIC :: pe_ffsl              ! component of pe_decomposed
  !
  ! Module variables
  !
  PUBLIC :: global_decomposition ! decomposition table for all PEs
  PUBLIC :: local_decomposition  ! decomposition info for this PE
  PUBLIC :: debug_parallel       ! Debug flag: -1= no debugging, 0,1=debugging
  PUBLIC :: debug_seriell        ! .true. to use old scheme to cycle longitudes
  PUBLIC :: any_col_1d           ! .true. if column model runs on any PE
  ! compatible with serial model version. 
  ! works only for nprocb==1
  !
  ! Module procedures
  !
  PUBLIC :: decompose              ! derive decomposition table
  PUBLIC :: print_decomposition    ! print decomposition table
  PUBLIC :: cleanup_decomposition  ! deallocate module variables
  !
  ! Information regarding the decomposition for 
  ! Flux-Form Semi-Lagrangian transport scheme

  TYPE pe_ffsl
    !
    SEQUENCE
    ! 
    INTEGER           :: nlat  ! number of latitudes on this pe                
    INTEGER           :: lats  ! southern  latitude  on this pe (S=1,N=nlatg)  
    INTEGER           :: latn  ! northern  latitudes on this pe                
    INTEGER           :: pe_s  ! processor id of neighbour in the south        
    INTEGER           :: pe_n  ! processor id of neighbour in the north        
    INTEGER           :: pe_x  ! processor to exchange hemispheres with        
    INTEGER           :: nlatx ! number of latitudes to exchange with pe_x     
  END TYPE pe_ffsl

  TYPE(pe_ffsl), PARAMETER :: default_ffsl = &
       pe_ffsl( &
       -1, & ! number of latitudes on this pe                       
       -1, & ! southern  latitude  on this pe (S=1,N=nlatg)  
       -1, & ! northern  latitudes on this pe                
       -1, & ! processor id of neighbour in the south        
       -1, & ! processor id of neighbour in the north        
       -1, & ! processor to exchange hemispheres with        
       -1)   ! number of latitudes to exchange with pe_x     

  ! Data type for local information of decomposition per PE 

  TYPE pe_decomposed
     !
     SEQUENCE
     !
     ! Information regarding the whole model domain and discretization

     INTEGER          :: nlon     ! number of longitudes
     INTEGER          :: nlat     ! number of latitudes
     INTEGER          :: npts = -1! number of global compressed grid points (e.g. land)
     LOGICAL, POINTER :: mask(:,:) => NULL() ! global mask
     INTEGER          :: nlev     ! number of levels
     INTEGER          :: nm       ! max. wavenumber used(triangular truncation)
     INTEGER ,POINTER :: nnp(:)   ! number of points on each m-column
     INTEGER ,POINTER :: nmp(:)   ! displacement of the first point of
                                  ! m-columns with respect to the first point 
                                  ! of the first m-column

     ! General information on the decomposition

     INTEGER          :: d_nprocs ! # of PEs for debugging/non-debugging domain
     INTEGER          :: spe      ! index # of first PE of domain
     INTEGER          :: epe      ! index # of last PE of domain
     INTEGER          :: nprocb   !#of PEs for dimension that counts longitudes
     INTEGER          :: nproca   !#of PEs for dimension that counts latitudes

     INTEGER ,POINTER :: mapmesh(:,:) ! indirection array mapping from a
                                      ! logical 2-d mesh to the processor index
                                      ! numbers

     ! local information

     INTEGER          :: pe        ! PE id 
     INTEGER          :: set_b     ! PE id in direction of longitudes
     INTEGER          :: set_a     ! PE id in direction of latitudes 
     LOGICAL          :: col_1d    ! 1d column model(s) on this pe

     ! Grid space decomposition 

     INTEGER          :: nglat     ! number of latitudes  on PE
     INTEGER          :: nglatmax  ! max. number of latitudes per PE
     INTEGER          :: nglon     ! number of longitudes on PE
     INTEGER          :: nglonmax  ! max. number of longitudes per PE
     INTEGER          :: nglh(2)   ! number of latitudes on each hemisphere
     INTEGER          :: glats(2)  ! start values of latitudes
     INTEGER          :: glate(2)  ! end values of latitudes
     INTEGER          :: glons(2)  ! start values of longitudes
     INTEGER          :: glone(2)  ! end values of longitudes
     INTEGER ,POINTER :: glat(:)   ! global latitude index N->S
     INTEGER ,POINTER :: glon(:)   ! offset to global longitude

     INTEGER          :: ngpts = -1! number of compressed grid points on PE
     INTEGER          :: gptss = -1! start value of grid indices
     INTEGER          :: gptse = -1! end value of grid indices

     ! Irregular grid

     INTEGER          :: ngpblks   ! number of rows
     INTEGER          :: nproma    ! number of columns
     INTEGER          :: npromz    ! number of columns in last row
     LOGICAL          :: lreg      ! flag for regular grid 

     ! Fourier space decomposition

     LOGICAL          :: lfused    ! true if this PE used in Fourier space
     INTEGER          :: nflat     ! number of latitudes on PE
     INTEGER          :: nflev     ! number of levels on PE
     INTEGER          :: nflevp1   ! number of levels+1 on PE
     INTEGER          :: flats(2)  ! start values of latitudes (on row)
     INTEGER          :: flate(2)  ! end values of latitudes (on row)
     INTEGER          :: flevs     ! start values of levels (on column)
     INTEGER          :: fleve     ! end values of levels (on column)

     ! Legendre space decomposition

     ! Row of PEs with same set_a
     INTEGER          :: nlm       ! number of local wave numbers m handled
     INTEGER ,POINTER :: lm(:)     ! actual local wave numbers m handled
     INTEGER          :: lnsp      ! number of complex spectral coefficients  
     INTEGER ,POINTER :: nlmp(:)   ! displacement of the first point of columns
     INTEGER ,POINTER :: nlnp(:)   ! number of points on each column
     INTEGER          :: nlnm0     ! number of coeff. with m=0 on this pe
     INTEGER ,POINTER :: intr(:)   ! index array used by transpose routine

     ! Column of PEs with same set_b
     INTEGER          :: nllev     ! number of levels
     INTEGER          :: nllevp1   ! number of levels+1
     INTEGER          :: llevs     ! start values of levels
     INTEGER          :: lleve     ! end values of levels

     ! Spectral space decomposition

     ! local PE
     INTEGER          :: snsp      ! number of spectral coefficients
     INTEGER          :: snsp2     ! 2*number of spectral coefficients
     INTEGER          :: ssps      ! first spectral coefficient
     INTEGER          :: sspe      ! last  spectral coefficient

     LOGICAL          :: lfirstc   ! true, if first global coeff (m=0,n=0) on
     ! this PE (for nudging)
     INTEGER          :: ifirstc   ! location of first global coeff on this PE
     INTEGER ,POINTER :: np1(:)    ! value of (n+1) for all coeffs of this PE 
     INTEGER ,POINTER :: np1i(:)   ! index of value of (n+1) for all coeffs
                                   ! of this PE 
     INTEGER ,POINTER :: mymsp(:)  ! value of m for all coeffs of this PE
     INTEGER          :: nns       ! number of different n-values for this PE
     INTEGER ,POINTER :: nindex(:) ! the nns elements contain the
     ! values of (n+1)

     INTEGER          :: nsm       ! number of wavenumbers per PE
     INTEGER ,POINTER :: sm (:)    ! actual local wave numbers handled
     INTEGER ,POINTER :: snnp(:)   ! number of n coeff. per wave number m
     INTEGER ,POINTER :: snn0(:)   ! first coeff. n for a given m
     INTEGER          :: nsnm0     ! number of coeffs with m=0 on this pe

     ! Flux-Form Semi-Lagrangian transport scheme decomposition

     TYPE (pe_ffsl)  :: ffsl
    
    ! Mark last element of derived type to allow proper usage later 
    ! in derived MPI datatypes 
    INTEGER            :: last          ! just dummy

  END TYPE pe_decomposed

! Fix for PGI compiler that doesn't handle named constants (default_ffsl) in
! type declaration with attribute PARAMETER
#if defined (__PGI)
  TYPE(pe_decomposed), SAVE :: default_decomposed = &
#else
  TYPE(pe_decomposed), PARAMETER :: default_decomposed = &
#endif
       pe_decomposed(    &
       0               , &
       0               , &
       -1              , &
       NULL()          , &
       0               , &
       0               , &
       NULL()          , &
       NULL()          , &
       0               , &
       0               , &
       0               , &
       0               , &
       0               , &
       NULL()          , &
       0               , &  
       0               , & 
       0               , &  
       .TRUE.          , &
       0               , & 
       0               , &
       0               , &
       0               , &
       (/ 0, 0 /)      , &
       (/ 0, 0 /)      , &
       (/ 0, 0 /)      , &
       (/ 0, 0 /)      , &
       (/ 0, 0 /)      , &
       NULL()          , &
       NULL()          , &
       -1              , &
       -1              , &
       -1              , &
       0               , &
       0               , &
       0               , & 
       .TRUE.          , &
       .TRUE.          , &
       0               , &
       0               , &
       0               , &
       (/ 0, 0 /)      , &
       (/ 0, 0 /)      , & 
       0               , &
       0               , &
       0               , &
       NULL()          , &
       0               , & 
       NULL()          , &
       NULL()          , &
       0               , &
       NULL()          , &
       0               , &
       0               , &
       0               , & 
       0               , & 
       0               , & 
       0               , & 
       0               , &
       0               , &
       .TRUE.          , &
       0               , &
       NULL()          , &
       NULL()          , &
       NULL()          , &
       0               , &
       NULL()          , &
       0               , &
       NULL()          , &
       NULL()          , &
       NULL()          , &
       0               , &
       default_ffsl    , &
       -1)
  
  ! Module variables

  TYPE (pe_decomposed), POINTER, SAVE :: global_decomposition(:)
  TYPE (pe_decomposed), SAVE          :: local_decomposition

   INTEGER                             :: debug_parallel = -1
  ! -1 no debugging
  ! 0 gather from PE0
  ! 1 gather from PE>0
  LOGICAL                             :: debug_seriell = .FALSE.
  ! .true. to use old scheme to cycle longitudes
  ! compatible with seriell model version.
  ! works only for nprocb==1
  LOGICAL                             :: any_col_1d = .FALSE.
  ! a column model is running on any of the PE's

  INTEGER                             :: iadvec
  ! the selected advection scheme

CONTAINS

  ! Module routines

  SUBROUTINE decompose (global_dc, nproma, nproca, nprocb, nlat, nlon, nlev, &
       nm, nn, nk, nadvec, mask, norot, debug, lfull_m, &
       lats_1d, lons_1d)
  !
  ! set decomposition table 
  !
  ! arguments
  !
  TYPE (pe_decomposed),INTENT(out) :: global_dc(0:)   ! decomposition table
  INTEGER             ,INTENT(in)  :: nproma          ! working dimension
  INTEGER             ,INTENT(in)  :: nproca, nprocb  ! no pe's in set A,B
  INTEGER             ,INTENT(in)  :: nlat, nlon, nlev! grid space size
  INTEGER             ,INTENT(in)  :: nk, nm, nn      ! truncation
  INTEGER             ,INTENT(in)  :: nadvec          ! advection scheme
  LOGICAL ,OPTIONAL   ,INTENT(in)  :: mask(:,:)       ! global mask
  LOGICAL ,OPTIONAL   ,INTENT(in)  :: norot           ! T:no lons rotation 
  INTEGER ,OPTIONAL   ,INTENT(in)  :: debug           ! run full model on PE0
  LOGICAL ,OPTIONAL   ,INTENT(in)  :: lfull_m         ! T: full columns
  INTEGER ,OPTIONAL   ,INTENT(in)  :: lats_1d         ! lats (column model)
  INTEGER ,OPTIONAL   ,INTENT(in)  :: lons_1d         ! lons (column model)
    !
    ! local variables
    !
    INTEGER :: pe, ii
    INTEGER :: nprocs
    LOGICAL :: nor
    LOGICAL :: lfullm
    !
    ! check consistency of arguments
    !
#ifndef STANDALONE
    IF (nlat/2 < nproca) THEN
       CALL finish ('decompose', &
       &            'Too many PEs selected for nproca - must be <= nlat/2.')
    END IF

    IF (nlat/nproca < 4) THEN
       CALL finish ('decompose', &
       &            'Too many PEs selected for nproca - nlat/nproca >= 4.')
    END IF
#endif

    IF (nlev < nprocb) THEN
       CALL finish ('decompose', &
       &            'Too many PEs selected for nprocb - must be < nlev.')
    END IF

    IF (nm+1 < nproca) THEN
       CALL finish ('decompose', &
       &            'Too many PEs selected for nproca - must be <= nm+1.')
    END IF

    IF (nk /= nm .OR. nn /= nm) THEN
       CALL finish ('decompose', &
       &            'Only triangular truncations supported in parallel mode.')
    END IF
    ! save advection scheme selector
    iadvec = nadvec
    !
    ! set local variables
    !
    nor    = .FALSE.; IF (PRESENT(norot )) nor = norot
    nprocs = nproca*nprocb
    IF (PRESENT(debug)) debug_parallel = debug
    lfullm = .FALSE.
    IF (PRESENT(lfull_m)) lfullm = lfull_m
    !
    ! set module variables
    !
    any_col_1d = PRESENT(lats_1d).AND.PRESENT(lons_1d)

    pe = 0
    IF (debug_parallel >= 0) THEN
       !
       ! set entries for full (debugging) model
       !
       ALLOCATE (global_dc(pe)%mapmesh(1:1,1:1))
       global_dc(pe)% d_nprocs = 1 ! # of PEs for debug domain
       global_dc(pe)% spe      = 1 ! index number of first PE for debug domain
       global_dc(pe)% epe      = 1 ! index number of last PE for debug domain
       !
       ! set entries for regular (partitioned) model
       !
       DO ii=pe+1,UBOUND(global_dc,1)
          ALLOCATE (global_dc(ii)%mapmesh(1:nprocb,1:nproca))
       END DO
       global_dc(pe+1:)% d_nprocs = p_nprocs - 1 ! normal domain
       global_dc(pe+1:)% spe      = 2            ! normal domain
       global_dc(pe+1:)% epe      = p_nprocs     ! normal domain
    ELSE
       !
       ! set entries for regular (partitioned) model
       !
       DO ii=pe,UBOUND(global_dc,1)
          ALLOCATE (global_dc(ii)%mapmesh(1:nprocb,1:nproca))
       END DO
       global_dc(pe:)% d_nprocs = p_nprocs
       global_dc(pe:)% spe      = 1
       global_dc(pe:)% epe      = p_nprocs
    END IF

    IF (debug_parallel >= 0) THEN
       !
       ! debug: PE 0 takes whole domain
       !
! Ha Ho-Hagemann {
       CALL set_decomposition (global_dc(pe:pe), 0, pe, nproca, nprocb, .TRUE., lfullm, mask=mask)
! Ha Ho-Hagemann }
!       CALL set_decomposition (global_dc(pe:pe), 0, pe, 1, 1, .TRUE., lfullm, mask=mask)
       ! from now on the coordinates of the debug PE in the logical
       ! mesh must be different from all other PEs' coordinates in order
       ! to avoid communication in the transposition routines
       !
       pe = pe + 1
       nprocs = nprocs + 1
    END IF
#if defined (__prism) && defined (use_comm_MPI1)
    CALL message('decompose','Could not check consistency of total PE number.')
#else
    !
    ! some more consistency checks
    !
    IF (nprocs /= SIZE(global_dc)) THEN
       WRITE(message_text,'(a,i0,i0)') &
            'nprocs, size(global_dc):', nprocs, SIZE(global_dc)
       CALL message('',message_text)
       CALL finish ('decompose', 'wrong size of global_dc.')
    END IF

    IF (nprocs /= p_nprocs) THEN
       CALL finish('decompose', &
            &            'Inconsistent total number of used PEs.')
    END IF
#endif
    !
    ! set entries for regular (partitioned) model
    !
    CALL set_decomposition (global_dc(pe:), nproma, pe, nproca, nprocb, nor, lfullm, &
                            mask=mask, lats_1d=lats_1d, lons_1d=lons_1d)
    !
    !
    !
    CALL p_set_communicator (nproca, nprocb, global_dc(pe)%mapmesh, & 
         debug_parallel)

  CONTAINS

    SUBROUTINE set_decomposition (dc, nproma, pe, nproca, nprocb, noro, lfullm, &
                                  mask, lats_1d, lons_1d)
    !
    ! set decomposition table for one model instance
    !
    TYPE (pe_decomposed) :: dc(0:)          ! decomposition table
    INTEGER ,INTENT(in)  :: nproma          ! working dimension
    INTEGER              :: pe              ! pe of dc(0)
    INTEGER              :: nproca, nprocb  ! no pe's in columns and rows
    LOGICAL              :: noro            ! T: no lons rotation
    LOGICAL              :: lfullm          ! T: full columns
    LOGICAL ,OPTIONAL    :: mask(:,:)       ! global mask
    INTEGER ,OPTIONAL    :: lats_1d         ! lats (column model)
    INTEGER ,OPTIONAL    :: lons_1d         ! lons (column model)

      INTEGER :: mepmash(nprocb,nproca)
      INTEGER :: i, p

      p = pe
      CALL decompose_global (dc, nlon, nlat, nlev, nm, nproca, nprocb, mask=mask)
      DO i = 0, UBOUND(dc,1)
         dc(i)%pe    = p                 ! local PE id 
         CALL mesh_map_init(-1, nprocb, nproca, pe, mepmash)
         dc(i)%mapmesh(:,:) = mepmash(:,:)

         CALL mesh_index(dc(i)%pe, nprocb, nproca, dc(i)%mapmesh, &
              dc(i)%set_b, dc(i)%set_a)

         CALL decompose_grid   (dc(i), nproma, noro, lats_1d=lats_1d, lons_1d=lons_1d)
         CALL decompose_level  (dc(i)) 
         CALL decompose_wavenumbers_row (dc(i))
         CALL decompose_specpoints_pe (dc(i), lfullm)

         p = p+1

      END DO

      IF (iadvec == 3) THEN
        DO i = 0, UBOUND(dc,1)
          CALL decompose_ffsl (dc, dc(i))
        END DO
      END IF

    END SUBROUTINE set_decomposition

  END SUBROUTINE decompose

  SUBROUTINE decompose_global (global_dc, nlon, nlat, nlev, nm, &
                               nproca, nprocb, mask)
  !--------------------------------------------------------
  ! set information on global domain in decomposition table
  !--------------------------------------------------------
  TYPE (pe_decomposed) ,INTENT(inout) :: global_dc(:) ! decomposition table
  INTEGER              ,INTENT(in)    :: nlon         ! # of longitudes
  INTEGER              ,INTENT(in)    :: nlat         ! # of latitudes
  INTEGER              ,INTENT(in)    :: nlev         ! # of levels
  INTEGER              ,INTENT(in)    :: nm           ! truncation
  INTEGER              ,INTENT(in)    :: nproca       ! # of PE's in N-S dir.
  INTEGER              ,INTENT(in)    :: nprocb       ! # of PE's in E-W dir.
  LOGICAL ,OPTIONAL    ,INTENT(in)    :: mask(:,:)    ! Global land mask
    !
    ! local scalars
    !
    INTEGER :: m, nn, nk, i
    !
    ! executable statements
    !
    global_dc% nlon   = nlon
    global_dc% nlat   = nlat
    global_dc% nlev   = nlev
    global_dc% nm     = nm
    global_dc% nproca = nproca
    global_dc% nprocb = nprocb

    IF (PRESENT(mask)) THEN
       global_dc% npts   = COUNT(mask)
       DO i=1,SIZE(global_dc)
          ALLOCATE(global_dc(i)%mask(nlon, nlat))
          global_dc(i)%mask = mask
       ENDDO
    ELSE
       global_dc% npts   = -1
    ENDIF

    nn = nm; nk = nm
    DO i=1,SIZE(global_dc,1)
       ALLOCATE (global_dc(i)% nnp(nm+1))
       ALLOCATE (global_dc(i)% nmp(nm+2))
       globaL_dc(i)%nnp=-999
       global_dc(i)%nmp=-999
    END DO
    DO m=0,nm
       global_dc(1)% nmp(1) = 0
       global_dc(1)% nnp(m+1) = MIN (nk-m,nn) + 1
       global_dc(1)% nmp(m+2) = global_dc(1)% nmp(m+1) + global_dc(1)% nnp(m+1)
    END DO
    DO i=2,SIZE(global_dc,1)
       global_dc(i)% nnp = global_dc(1)% nnp
       global_dc(i)% nmp = global_dc(1)% nmp
    END DO

  END SUBROUTINE decompose_global

  !============================================================================

  SUBROUTINE decompose_grid (pe_dc, nproma, norot, lats_1d, lons_1d) 
  !
  ! sets information on local gridpoint space in decomposition table entry
  !
  TYPE (pe_decomposed) ,INTENT(inout) :: pe_dc ! decomposition table entry
  INTEGER              ,INTENT(in)    :: nproma  ! working dimension
  LOGICAL              ,INTENT(in)    :: norot   ! rotate domain in south. hem.
  INTEGER  ,OPTIONAL   ,INTENT(in)    :: lats_1d ! latitude for column model
  INTEGER  ,OPTIONAL   ,INTENT(in)    :: lons_1d ! longitude for column model
    !
    ! local variables
    !
    INTEGER :: ngl, i, inp, inprest
    INTEGER :: nptrlat(pe_dc% nproca+1), nptrlon(pe_dc% nprocb+1), nptr(pe_dc%nproca*pe_dc%nprocb+1)
    INTEGER :: set_a, set_b, nlon, nlat, npts, nproca, nprocb
    !
    ! executable statements
    !
    IF(PRESENT(lats_1d).AND.PRESENT(lons_1d)) THEN
      !
      ! column model settings
      !
      pe_dc% col_1d = .TRUE.
      !
      ! latitudes
      !
      pe_dc% nglat    = 1
      pe_dc% nglh     = (/1,0/)
      pe_dc% glats(1) = lats_1d
      pe_dc% glate(1) = lats_1d
      pe_dc% glats(2) = 1
      pe_dc% glate(2) = 0
      ALLOCATE (pe_dc% glat(1))
      pe_dc% glat     = lats_1d
      !
      ! longitudes
      !
      pe_dc% nglon    = 1
      pe_dc% glons(1) = lons_1d
      pe_dc% glone(1) = lons_1d
      pe_dc% glons(2) = 1
      pe_dc% glone(2) = 0
      ALLOCATE (pe_dc% glon(1))
      pe_dc% glon     = lons_1d - 1

      pe_dc% nproma = 1
      pe_dc% ngpblks = 1
      pe_dc% npromz = 0

    ELSE IF(pe_dc%npts > 0) THEN
      pe_dc%col_1d = .FALSE.
      pe_dc%lreg = .FALSE.

      set_a  = pe_dc% set_a
      set_b  = pe_dc% set_b
      nproca = pe_dc% nproca
      nprocb = pe_dc% nprocb
      npts   = pe_dc% npts

      nptr(:) = -999
      inp = npts / (nproca*nprocb)
      inprest = npts - inp*nproca*nprocb
      nptr(1) = 1
      DO i=1,nproca*nprocb
         IF (i <= inprest) THEN
            nptr(i+1) = nptr(i) + inp + 1
         ELSE
            nptr(i+1) = nptr(i) + inp
         ENDIF
      ENDDO

      i = (set_a - 1)*nprocb + set_b
      pe_dc%ngpts = nptr(i+1) - nptr(i)
      pe_dc%gptss = nptr(i)
      pe_dc%gptse = nptr(i+1) - 1

      IF (nproma > 0) THEN
         pe_dc%nproma = MIN(nproma, pe_dc%ngpts)
      ELSE
         pe_dc%nproma = pe_dc%ngpts
      ENDIF
      pe_dc% ngpblks = ( pe_dc% ngpts - 1 ) / pe_dc% nproma + 1
      pe_dc% npromz = pe_dc% ngpts  - (pe_dc% ngpblks - 1) * pe_dc% nproma

      ! local copies of decomposition table entries
      !
      nlon   = pe_dc% nlon
      nlat   = pe_dc% nlat

      pe_dc% nglat    =  0
      pe_dc% nglh     =  0
      pe_dc% glats(1) = 0
      pe_dc% glate(1) = 0
      pe_dc% glats(2) = 0
      pe_dc% glate(2) = 0
      pe_dc% nglon    = 0
      pe_dc% glons(1) = 0
      pe_dc% glone(1) = 0
      NULLIFY (pe_dc% glat)
      NULLIFY (pe_dc% glon)

    ELSE
      !
      ! full model settings
      !
      pe_dc% col_1d = .FALSE.
      !
      ! local copies of decomposition table entries
      !
      set_a  = pe_dc% set_a
      set_b  = pe_dc% set_b
      nlon   = pe_dc% nlon
      nlat   = pe_dc% nlat
      nproca = pe_dc% nproca
      nprocb = pe_dc% nprocb
      !
      ! first distribute latitudes
      !
      ngl = nlat/2

      nptrlat(:) = -999 
      pe_dc% nglatmax = 0

      inp = ngl/nproca
      inprest = ngl-inp*nproca
      nptrlat(1) = 1
      DO i = 1, nproca
         IF (i <= inprest) THEN
            nptrlat(i+1) = nptrlat(i)+inp+1
         ELSE
            nptrlat(i+1) = nptrlat(i)+inp   
         END IF
         pe_dc% nglatmax =MAX(pe_dc% nglatmax,2*(nptrlat(i+1)-nptrlat(i)))
      END DO
      !
      ! now distribute longitudes
      !
      nptrlon(:) = -999 
      pe_dc% nglonmax = 0

      inp = nlon/nprocb
      inprest = nlon-inp*nprocb
      nptrlon(1) = 1
      DO i = 1, nprocb
         IF (i <= inprest) THEN
            nptrlon(i+1) = nptrlon(i)+inp+1
         ELSE
            nptrlon(i+1) = nptrlon(i)+inp   
         END IF
         pe_dc% nglonmax =MAX(pe_dc% nglonmax,nptrlon(i+1)-nptrlon(i))
      END DO

      ! define parts per pe
      !
      ! latitudes
      !
      pe_dc% nglat    =  2*(nptrlat(set_a+1)-nptrlat(set_a))
      pe_dc% nglh     =  pe_dc% nglat/2

      pe_dc% glats(1) = nptrlat(set_a)
      pe_dc% glate(1) = nptrlat(set_a+1)-1

      pe_dc% glats(2) = nlat-nptrlat(set_a+1)+2 
      pe_dc% glate(2) = nlat-nptrlat(set_a)+1

      ALLOCATE (pe_dc% glat( pe_dc% nglat))
      pe_dc% glat(1:pe_dc% nglat/2) = (/(i,i=pe_dc% glats(1),pe_dc% glate(1))/)
      pe_dc% glat(pe_dc% nglat/2+1:)= (/(i,i=pe_dc% glats(2),pe_dc% glate(2))/)
      !
      ! longitudes
      !
      pe_dc% nglon    = nptrlon(set_b+1)-nptrlon(set_b)
      pe_dc% glons(1) = nptrlon(set_b)
      pe_dc% glone(1) = nptrlon(set_b+1)-1
      !
      ! rotate longitudes in southern area
      !
      IF (norot) THEN
         pe_dc% glons(2) = pe_dc% glons(1)
         pe_dc% glone(2) = pe_dc% glone(1)
      ELSE
         pe_dc% glons(2) = MOD(nptrlon(set_b)-1+nlon/2, nlon)+1
         pe_dc% glone(2) = MOD(nptrlon(set_b+1)-2+nlon/2, nlon)+1
      END IF

      ALLOCATE (pe_dc% glon( pe_dc% nglat))
      pe_dc% glon(1:pe_dc% nglat/2)  = pe_dc% glons(1)-1
      pe_dc% glon(pe_dc% nglat/2+1:) = pe_dc% glons(2)-1

    ! setting for irregular grid

    IF ( nproma > 0 ) THEN
      pe_dc% lreg   = .FALSE.
      pe_dc% nproma = MIN (nproma, pe_dc% nglon * pe_dc% nglat)
    ELSE
      pe_dc% lreg   = .TRUE.
      pe_dc% nproma = pe_dc% nglon
    ENDIF

    pe_dc% ngpblks = ( pe_dc% nglon * pe_dc% nglat - 1 ) / pe_dc% nproma + 1

    pe_dc% npromz = pe_dc% nglon * pe_dc% nglat - (pe_dc% ngpblks - 1) * pe_dc% nproma

    ENDIF

  END SUBROUTINE decompose_grid

  SUBROUTINE decompose_level (pe_dc)
  !
  ! determine number of local levels for fourier and legendre calculations
  ! this is based on the supplied nlev and nprocb
  !
  TYPE (pe_decomposed), INTENT(inout) :: pe_dc

    INTEGER :: set_b, nlev, nprocb
    INTEGER :: inp, inprest, jb

    INTEGER :: ll(pe_dc% nprocb+1), nptrll(pe_dc% nprocb+1)
    !
    ! copy table entries to local variables
    !
    set_b  = pe_dc% set_b
    nlev   = pe_dc% nlev
    nprocb = pe_dc% nprocb
    IF (pe_dc% col_1d) THEN
      !
      ! settings for column model
      !
      pe_dc% nflat   = 0
      pe_dc% flats   = 1
      pe_dc% flate   = 0
      pe_dc% nflevp1 = 0
      pe_dc% flevs   = 1
      pe_dc% fleve   = 0
      pe_dc% lfused  = .FALSE.
      pe_dc% nflev   = 0
    ELSE
      !
      ! settings for full model
      !
      ! latitudes are the same as in grid space
      !
      pe_dc% nflat = pe_dc% nglat
      pe_dc% flats = pe_dc% glats
      pe_dc% flate = pe_dc% glate
      !
      ! distribute levels
      !
      inp = (nlev+1)/nprocb
      inprest = (nlev+1)-inp*nprocb
      !
      ! make sure highest level is on a PE which has one of the subsets with an
      ! extra level => improves load-balance
      !
      nptrll(nprocb+1)=nlev+1+1
      DO jb = nprocb,1,-1
         IF ((nprocb - jb + 1) <= inprest) THEN 
            ll(jb) = inp+1
         ELSE
            ll(jb) = inp
         END IF
         nptrll(jb) = nptrll(jb+1)-ll(jb)
      END DO
      pe_dc% nflevp1 = ll(set_b)
      pe_dc% flevs   = nptrll(set_b)
      pe_dc% fleve   = nptrll(set_b+1)-1
      pe_dc% lfused  = pe_dc%nflevp1 > 0

      IF ( pe_dc%fleve > nlev ) THEN
         pe_dc%nflev = pe_dc%nflevp1 - 1
      ELSE
         pe_dc%nflev = pe_dc%nflevp1
      END IF
    ENDIF
  END SUBROUTINE decompose_level

  SUBROUTINE decompose_wavenumbers_row (pe_dc)
  TYPE (pe_decomposed), INTENT(inout) :: pe_dc
  !
  ! decompose wavenumbers along latitudinal direction (nproca)
  !
    INTEGER :: jm, ik, il, ind

    INTEGER :: set_a, nm, nproca

    INTEGER :: nprocm(0:pe_dc% nm)   ! 'set_a' for a certain spectral wave
    INTEGER :: nspec (pe_dc% nproca) ! complex spectral coeff. per row (set_a)
    INTEGER :: numpp (pe_dc% nproca) ! spectral waves per processor row (set_a)
    INTEGER :: myms  (pe_dc% nm+1)   ! actual wave numbers handled
    !
    ! copy decomposition table entries to local variables
    !
    set_a  = pe_dc% set_a
    nm     = pe_dc% nm
    nproca = pe_dc% nproca
    !
    ! levels same as in Fourier space
    !
    pe_dc% nllev   = pe_dc% nflev
    pe_dc% nllevp1 = pe_dc% nflevp1
    pe_dc% llevs   = pe_dc% flevs
    pe_dc% lleve   = pe_dc% fleve
    !
    ! distribute spectral waves
    !
    nspec(:) = 0
    numpp(:) = 0
    myms(:)  = -1
    il       = 1

    IF (pe_dc% col_1d) THEN
      !
      ! column model settings
      !
      nm = -1
    ELSE
      !
      ! full model settings
      !
      ind = 1
      ik  = 0
      DO jm = 0, nm
         ik = ik + ind
         IF (ik > nproca) THEN
            ik = nproca
            ind = -1
         ELSE IF (ik < 1) THEN
            ik = 1
            ind = 1
         END IF
         nprocm(jm) = ik
         nspec(ik) = nspec(ik)+nm-jm+1
         numpp(ik) = numpp(ik)+1
         IF (ik == set_a) THEN
            myms(il) = jm
            il = il+1
         END IF
      END DO
    ENDIF

    ALLOCATE(pe_dc% lm   (numpp(set_a)  ))
    ALLOCATE(pe_dc% nlmp (numpp(set_a)+1))
    ALLOCATE(pe_dc% nlnp (numpp(set_a)  ))
    ALLOCATE(pe_dc% intr (numpp(set_a)*2))

    pe_dc% nlm     = numpp(set_a)
    pe_dc% lnsp    = nspec(set_a)
    pe_dc% lm      = myms (1:il-1)
    pe_dc% nlnp    = nm - pe_dc% lm + 1
    pe_dc% nlnm0   = 0; IF (myms (1)== 0) pe_dc% nlnm0 = pe_dc% nlnp(1)
    pe_dc% nlmp(1) = 0
    DO jm = 1, pe_dc% nlm
       pe_dc% nlmp(  jm+1) = pe_dc% nlmp(jm) + pe_dc% nlnp(jm)
       pe_dc% intr(2*jm-1) = 2 * pe_dc% lm(jm) + 1
       pe_dc% intr(2*jm  ) = 2 * pe_dc% lm(jm) + 2
    END DO

  END SUBROUTINE decompose_wavenumbers_row

  SUBROUTINE decompose_specpoints_pe (pe_dc, lfullm)
  LOGICAL ,INTENT(in) :: lfullm

    ! Decompose wavenumbers additionally along longitudinal direction (nprocb),
    ! thus decompose wavenumbers among individual PEs.
    ! Each PE receives full or partial m-columns, depending on strategy 
    ! chosen. The case with partial m-columns will generally be better load-
    ! balanced.


    TYPE (pe_decomposed), INTENT(inout) :: pe_dc

    INTEGER :: nspec, set_b, nump, myml(pe_dc%nlm)

    INTEGER :: insef, irestf, jb, jmloc, icolset, iave, im, inm, iii, jn, i
    INTEGER :: ispec, il, jm, imp, inp, is, ic, inns, nsm, ifirstc
    INTEGER :: nm, nprocb

    INTEGER :: nptrsv (pe_dc% nprocb+1) ! first spectral wave column (PE)
    INTEGER :: nptrsvf(pe_dc% nprocb+1) ! full spectral wave m-columns (PE)
    INTEGER :: nptrmf (pe_dc% nprocb+1) ! distribution of m-columns among PE
    ! columns (used for semi-impl. calculations in full m-columns case)
    INTEGER :: nspstaf(0:pe_dc% nm)     ! m-column starts here (used for 
    ! semi-impl. full m-columns case)

    INTEGER :: inumsvf(pe_dc% nprocb+1)

    INTEGER :: nns               ! number of different n-values for this PE
    INTEGER :: nindex(pe_dc%nm+1)! the first nns elements contain the values 
    ! of (n+1);
    ! for non triangular truncation: pe_dc%nkp1
    LOGICAL :: lnp1(pe_dc%nm+1)  ! if false: this value of (n+1) has not 
    ! occurred yet
    INTEGER :: itmp(pe_dc%nm+1)

    INTEGER :: np1(pe_dc%lnsp)   ! np1 holds the values of (n+1) of all
    ! coeffs on this PE
    INTEGER :: mymsp(pe_dc%lnsp) ! mymsp holds the values of m of all
    ! coeffs on this PE

    INTEGER :: myms(pe_dc%nlm)   ! myms holds the values of the wavenumbers
    ! on this PE 
    INTEGER :: myns(pe_dc%nlm)   ! myns holds the number of coeffs of
    ! m-columns (full or partial) on this PE
    INTEGER :: snn0(pe_dc%nlm)   ! snn0 holds the offset of n wavenumbers
    ! for each column (0 for full)

    LOGICAL :: ljm, lfirstc

    ! executable statements

    ! short names for processor row variables

    nspec  = pe_dc% lnsp   ! number of spectral coefficients in processor row
    nump   = pe_dc% nlm    ! number of wavenumbers in processor row
    myml   = pe_dc% lm     ! wave numbers handled in processor row
    nm     = pe_dc% nm     ! max. wavenumber used(triangular truncation)
    nprocb = pe_dc% nprocb ! number of PEs in set
    set_b  = pe_dc% set_b  ! index of PE in set

    ! Partitioning of spectral coefficients in semi-implicit calculations.
    ! Be careful : nspec varies between processor rows, so nptrsv() 
    ! differs on different processor rows.

    insef = nspec/nprocb
    irestf = nspec-insef*nprocb
    nptrsv(1) = 1
    DO jb = 2, nprocb+1
       IF(jb-1 <= irestf) THEN
          nptrsv(jb) = nptrsv(jb-1)+insef+1
       ELSE
          nptrsv(jb) = nptrsv(jb-1)+insef
       END IF
    END DO

    ! partitioning of spectral coefficients in semi-implicit calculations
    ! for the case where complete m-columns are required.

    ! original idea

    nptrmf(:) = 1
    inumsvf(:) = 0
    icolset = MIN(nprocb,nump)
    nptrmf(1) = 1
    nptrmf(icolset+1:nprocb+1) = nump+1
    ispec = nspec
    IF (nump>0) THEN
      iave = (ispec-1)/icolset+1
      DO jmloc = nump, 1, -1
         im = myml(jmloc)
         inm = nm-im+1
         IF (inumsvf(icolset) < iave) THEN
            inumsvf(icolset) = inumsvf(icolset)+inm
         ELSE
            nptrmf(icolset) = jmloc+1
            ispec = ispec-inumsvf(icolset)
            icolset = icolset-1
            IF (icolset == 0) THEN
               CALL finish ('decompose_wavenumbers_pe', &
                    &                'error in decomposition, icolset = 0!')
            END IF
            iave = (ispec-1)/icolset+1
            inumsvf(icolset) = inumsvf(icolset)+inm
         END IF
      END DO
    ENDIF
    nptrsvf(1) = 1
    DO jb = 2, nprocb+1
       nptrsvf(jb) = nptrsvf(jb-1)+inumsvf(jb-1)
    END DO

    nspstaf(:) = -999
    iii = 1
    !    DO jmloc = nptrmf(set_b), nptrmf(set_b+1)-1
    !      nspstaf(myml(jmloc)) = iii
    !      iii = iii+nm-myms(jmloc)+1
    !    END DO

    IF (lfullm) THEN
       pe_dc%ssps = nptrsvf(set_b)
       pe_dc%sspe = nptrsvf(set_b+1)-1
       pe_dc%snsp = nptrsvf(set_b+1)-nptrsvf(set_b)
       pe_dc%snsp2= 2*pe_dc%snsp
    ELSE
       pe_dc%ssps = nptrsv(set_b)
       pe_dc%sspe = nptrsv(set_b+1)-1
       pe_dc%snsp = nptrsv(set_b+1)-nptrsv(set_b)
       pe_dc%snsp2= 2*pe_dc%snsp
    END IF

    ! il counts spectral coeffs per PE
    ! nsm counts number of wavenumbers (full or partial) per PE

    il  = 0
    nsm = 0
    ljm = .FALSE.
    lnp1(:) = .FALSE.
    nindex(:) = -999
    inns = 0
    myms(:) = -999
    myns(:) = 0
    ifirstc = -999
    lfirstc = .FALSE.

    ! loop over all coeffs of a processor row and check whether they belong to
    ! this PE

    DO jm = 1,nump
       imp = pe_dc%nlmp(jm)
       inp = pe_dc%nlnp(jm)
       DO jn = 1,inp
          is = imp + jn
          ic = myml(jm) +jn

          IF ((is >= pe_dc%ssps) .AND. (is <= pe_dc%sspe)) THEN
             IF (.NOT.ljm) snn0 (nsm+1) = jn-1
             ljm = .TRUE.
             il  = il + 1

             np1(il) = ic
             mymsp(il) = myml(jm)

             myns(nsm+1) = myns(nsm+1) + 1

             IF (.NOT. lnp1(ic)) THEN
                inns = inns + 1
                nindex(inns) = ic
                lnp1(ic) = .TRUE.
             END IF

             ! first global coeff. on this PE? (for nudging)
             IF ((mymsp(il) == 0) .AND. (np1(il) == 1)) THEN
                lfirstc = .TRUE.
                ifirstc = il
             END IF

          END IF
       END DO
       IF (ljm) THEN
          nsm = nsm + 1
          myms(nsm) = myml(jm)
       END IF
       ljm = .FALSE.
    END DO

    nns = inns

    IF (il /= pe_dc%snsp) THEN
       CALL finish('decompose', &
            &            'Error in computing number of spectral coeffs')
    END IF

    ALLOCATE (pe_dc%np1(il))
    ALLOCATE (pe_dc%np1i(il))
    ALLOCATE (pe_dc%mymsp(il))
    pe_dc%np1 = np1(1:il)

    DO i=1,nns
      itmp(nindex(i))=i
    ENDDO
    DO i=1,il
      pe_dc%np1i(i)=itmp(np1(i))
    ENDDO

    pe_dc%mymsp = mymsp(1:il)

    pe_dc%nns = nns
    ALLOCATE (pe_dc%nindex(pe_dc%nns))
    pe_dc%nindex(:) = nindex(1:nns)

    pe_dc%lfirstc = lfirstc
    pe_dc%ifirstc = ifirstc

    pe_dc%nsm = nsm
    ALLOCATE (pe_dc%sm (pe_dc%nsm))
    ALLOCATE (pe_dc%snnp (pe_dc%nsm))
    pe_dc%sm = myms(1:nsm)
    pe_dc%snnp = myns(1:nsm)

    ALLOCATE (pe_dc%snn0 (pe_dc%nsm))
    pe_dc%snn0 = snn0 (1:nsm)

    pe_dc%nsnm0 = 0
    DO i=1,nsm
       IF (pe_dc%sm(i) == 0) THEN
          pe_dc%nsnm0 = pe_dc%snnp(i)
          EXIT
       END IF
    END DO

  END SUBROUTINE decompose_specpoints_pe

  SUBROUTINE print_decomposition (dc)

    TYPE (pe_decomposed), INTENT(in) :: dc

    CALL message('', '')
    CALL message('', 'Detailed ')
!OSBUTLS    CALL message('', separator)
    IF (.NOT. ASSOCIATED(dc%lm)) THEN
       CALL message('print_decomposition', &
            'decomposition not done, cannot print information ...')
    END IF

    WRITE (message_text,'(a,i5)') ' PE    : ', dc%pe
    CALL message('', message_text)
    WRITE (message_text,'(a,i5)') ' Processor row    (Set A) : ', dc%set_a
    CALL message('', message_text)
    WRITE (message_text,'(a,i5)') ' Processor column (Set B) : ', dc%set_b
    CALL message('', message_text)
    WRITE (message_text,'(a,l5)') ' Column model running     : ', dc%col_1d
    CALL message('', message_text)

    CALL message('','')

    WRITE (message_text,'(a)') ' mapmesh : '
    CALL message('', message_text)
    WRITE (message_text,'(12i5)') dc%mapmesh
    CALL message('', message_text)

    CALL message('','')

    WRITE (message_text,'(a,i5)') ' spe : ', dc%spe
    CALL message('', message_text)
    WRITE (message_text,'(a,i5)') ' epe : ', dc%epe
    CALL message('', message_text)
    WRITE (message_text,'(a,i5)') ' d_nprocs : ', dc%d_nprocs
    CALL message('', message_text)

    CALL message('','')

    WRITE (message_text,'(a,i5)') ' nlon  : ', dc%nlon
    CALL message('', message_text)
    WRITE (message_text,'(a,i5)') ' nlat  : ', dc%nlat
    CALL message('', message_text)
    WRITE (message_text,'(a,i5)') ' nlev  : ', dc%nlev
    CALL message('', message_text)
    WRITE (message_text,'(a,i5)') ' nm    : ', dc%nm
    CALL message('', message_text)
    CALL print_field(' nmp   : ', dc%nmp)
    CALL print_field(' nnp   : ', dc%nnp)
    WRITE (message_text,'(a,i5)') ' nproca: ', dc%nproca
    CALL message('', message_text)
    WRITE (message_text,'(a,i5)') ' nprocb: ', dc%nprocb
    CALL message('', message_text)
    WRITE (message_text,'(a,i5)') ' npts:   ', dc%npts
    CALL message('', message_text)

    CALL message('','')

    WRITE (message_text,'(i5,a)') dc%nglat, ' latitudes  in grid space'
    CALL message('', message_text)
    WRITE (message_text,'(i5,a)') dc%nglon, ' longitudes in grid space'
    CALL message('', message_text)
    WRITE (message_text,'(2i5,a)')dc%nglh,  ' latitudes on N/S hemisphere'
    CALL message('', message_text)
    WRITE (message_text,'(i5,a)') dc%nglat*dc%nglon, ' grid points (total)'
    CALL message('', message_text)
    WRITE (message_text,'(4(a,i5))') ' glatse: ', dc%glats(1), ' -', dc%glate(1), &
         ',   ', dc%glats(2), ' -', dc%glate(2)
    CALL message('', message_text)
    WRITE (message_text,'(4(a,i5))') ' glonse: ', dc%glons(1), ' -', dc%glone(1), &
         ',   ', dc%glons(2), ' -', dc%glone(2)
    CALL message('', message_text)
    IF (ASSOCIATED(dc%glat)) CALL print_field(' glat    : ', dc% glat)
    IF (ASSOCIATED(dc%glon)) CALL print_field(' glon    : ', dc% glon)
    IF (dc%npts /= -1) THEN
       WRITE (message_text,'(2(a,i5))') ' gptsse: ', dc%gptss, ' -', dc%gptse
       CALL message('', message_text)
    ENDIF

    CALL message('','')

    WRITE (message_text,'(l4,a)') dc%lfused,  ' processor used in Fourier space'
    CALL message('', message_text)
    WRITE (message_text,'(i5,a)') dc%nflat,   ' latitudes  in Fourier space'
    CALL message('', message_text)
    WRITE (message_text,'(i5,a)') dc%nflevp1, ' levels+1 in Fourier space'
    CALL message('', message_text)
    WRITE (message_text,'(i5,a)') dc%nflev,   ' levels in Fourier space'
    CALL message('', message_text)
    WRITE (message_text,'(4(a,i5))') ' flat  : ', dc%flats(1), ' -', dc%flate(1), &
         ',   ', dc%flats(2), ' -', dc%flate(2)
    CALL message('', message_text)
    WRITE (message_text,'(2(a,i5))') ' flev  : ', dc%flevs,    ' -', dc%fleve
    CALL message('', message_text)

    CALL message('','')

    WRITE (message_text,'(i5,a)') dc%nlm,     ' wave numbers in Legendre space'
    CALL message('', message_text)
    CALL print_field(' lm    : ', dc%lm(:dc%nlm))
    CALL print_field(' nlnp  : ', dc%nlnp(:dc%nlm))
    CALL print_field(' nlmp  : ', dc%nlmp(:dc%nlm+1))
    CALL message('', message_text)
    WRITE (message_text,'(i5,a)') dc%nlnm0,   ' coefficients for m=0'
    CALL message('', message_text)
    WRITE (message_text,'(i5,a)') dc%lnsp,    ' spectral coefficients'
    CALL message('', message_text)
    WRITE (message_text,'(i5,a)') dc%nllevp1, ' levels+1 in Legendre space'
    CALL message('', message_text)
    WRITE (message_text,'(i5,a)') dc%nllev,   ' levels in Legendre space'
    CALL message('', message_text)
    WRITE (message_text,'(2(a,i5))')          ' llev  : ', dc%llevs, ' -', dc%lleve
    CALL message('', message_text)

    CALL message('','')

    WRITE (message_text,'(i5,a)') dc%snsp,   ' coefficients'
    CALL message('', message_text)
    WRITE (message_text,'(i5,a)') dc%snsp2,  ' coefficients times 2'
    CALL message('', message_text)
    WRITE (message_text,'(i5,a)') dc%ssps,   ' first spectral coefficient'
    CALL message('', message_text)
    WRITE (message_text,'(i5,a)') dc%sspe,   ' last spectral coefficient'
    CALL message('', message_text)

    WRITE (message_text,'(l4,a)') dc%lfirstc, ' first global coefficient on this PE'
    CALL message('', message_text)
    WRITE (message_text,'(i5,a)') dc%ifirstc, ' local index of first global coefficient'
    CALL message('', message_text)

    CALL print_field(' np1   : ', dc%np1) 
    CALL print_field(' mymsp : ', dc%mymsp)
    WRITE (message_text,'(i5,a)') dc%nns ,    ' number of different n-values' 
    CALL message('', message_text)
    CALL print_field(' nindex: ', dc%nindex)

    WRITE (message_text,'(i5,a)') dc%nsm ,    ' number of m wave numbers.'
    CALL message('', message_text)
    CALL print_field(' sm    : ', dc%sm) 
    CALL print_field(' snnp  : ', dc%snnp)
    CALL print_field(' snn0  : ', dc%snn0)
    WRITE (message_text,'(i5,a)') dc%nsnm0 ,  ' coefficients for m=0'
    CALL message('', message_text)

    IF (iadvec == 3) THEN
       CALL print_ffsl_dec ((/dc/))
    ENDIF

    CONTAINS

      SUBROUTINE print_field(title, field)
        
        CHARACTER(len=*), INTENT(in) :: title
        INTEGER,          INTENT(in) :: field(:)
        
        INTEGER :: i, isize
        
        CHARACTER(len=16) :: cformat

        isize = SIZE(field)

        IF (isize <= 12) THEN       
          WRITE(cformat,'(a,i0,a)') '(a,', isize, 'i5)'
          WRITE (message_text,TRIM(cformat)) title, field(1:isize)
          CALL message ('',message_text)
        ELSE
          WRITE (message_text,'(a,12i5)') title, field(1:12)
          DO i = 13, isize, 12
            IF (isize <= i+11) THEN
              WRITE(cformat,'(a,i0,a)') '(9x,', isize, 'i5)'
              WRITE (message_text,TRIM(cformat)) field(i:isize)
              CALL message ('',message_text)
            ELSE
              WRITE (message_text,'(9x,12i5)') field(i:i+11)
              CALL message ('',message_text)
            ENDIF
          ENDDO
        ENDIF
        
      END SUBROUTINE print_field
      
  END SUBROUTINE print_decomposition

  SUBROUTINE mesh_map_init(option, isize, jsize, p, meshmap)

    ! all isize*jsize PEs are mapped to a 2-D logical mesh
    ! starting with PE with id p

    INTEGER :: option, isize, jsize, p
    INTEGER :: meshmap(1:isize,1:jsize)
    INTEGER :: i, j

    IF (option == 1) THEN
       ! row major ordering
       DO j = 1, jsize
          DO i = 1, isize
             meshmap(i,j) = (i-1) + (j-1)*isize + p
          END DO
       END DO
    ELSE IF (option == -1) THEN
       !column major ordering
       DO j = 1, jsize
          DO i = 1, isize
             meshmap(i,j) = (j-1) + (i-1)*jsize + p
          END DO
       END DO
    END IF

  END SUBROUTINE mesh_map_init

  SUBROUTINE mesh_index(p, isize, jsize, meshmap, idex, jdex)

    ! returns coordinates (idex,jdex) of PE with id p in logical mesh

    INTEGER :: p, isize, jsize, idex, jdex
    INTEGER :: meshmap(1:isize,1:jsize)
    INTEGER :: i, j
    LOGICAL :: lexit

    idex  = -1
    jdex  = -1
    lexit = .FALSE.

    DO j = 1, jsize
       DO i = 1, isize
          IF (meshmap(i,j) == p) THEN
             ! coordinates start at 1
             idex  = i
             jdex  = j
             lexit = .TRUE.
             EXIT
          END IF
       END DO
       IF (lexit) EXIT
    END DO
    ! check for successful completion of search
    IF ((idex == -1) .OR. (jdex == -1)) THEN
       CALL finish('mesh_index', &
            &            'Unable to find processor in meshmap array')
    END IF

  END SUBROUTINE mesh_index
!==============================================================================
  SUBROUTINE decompose_ffsl (dc, pe_dc)
  TYPE(pe_decomposed) ,INTENT (in)    :: dc (:) ! global decomposition array
  TYPE(pe_decomposed) ,INTENT (inout) :: pe_dc  ! local  decomposition element 
  !
  ! decomposition for Flux-Form Semi-Lagrangian transport scheme
  !
    INTEGER :: set_a, set_b, nproca, nprocb
    INTEGER :: i
    set_a  = pe_dc% set_a
    set_b  = pe_dc% set_b
    nproca = pe_dc% nproca
    nprocb = pe_dc% nprocb
    !
    ! currently only works for non rotated southern hemisphere
    ! levels may be distributed later
    !
    IF (pe_dc% glons(1) /= pe_dc% glons(2)) &
                   CALL finish ('decompose_ffsl','rotated grid')

    pe_dc% ffsl% pe_x = -1
    pe_dc% ffsl% pe_n = -1
    pe_dc% ffsl% pe_s = -1        

    IF(nproca==1) THEN
      !
      ! settings for single processor or column model mode
      !
      pe_dc% ffsl% pe_x  = pe_dc% pe
      pe_dc% ffsl% nlat  = pe_dc% nglat
      pe_dc% ffsl% nlatx = pe_dc% nglat/2
      pe_dc% ffsl% lats  = pe_dc% glats(1) + pe_dc% nglat - 1
      pe_dc% ffsl% latn  = pe_dc% glats(1)
    ELSE
      !
      ! settings for parallel mode
      !
      IF(MOD(set_a,2)==1) THEN        ! keep northern hemisphere on this pe
        IF(set_a == nproca) THEN
          pe_dc% ffsl% pe_x = pe_dc% pe
          pe_dc% ffsl% pe_s = pe_dc% mapmesh (set_b,set_a-1)
        ELSE IF(set_a == nproca-1) THEN
          pe_dc% ffsl% pe_x = pe_dc% mapmesh (set_b,set_a+1)
          pe_dc% ffsl% pe_s = pe_dc% mapmesh (set_b,set_a+1)
        ELSE
          pe_dc% ffsl% pe_x = pe_dc% mapmesh (set_b,set_a+1)
          pe_dc% ffsl% pe_s = pe_dc% mapmesh (set_b,set_a+2)
        ENDIF
        IF(set_a /= 1) THEN
          pe_dc% ffsl% pe_n = pe_dc% mapmesh (set_b,set_a-2)
        ENDIF
        pe_dc% ffsl% latn = pe_dc% glats (1)
        i = indx (pe_dc% ffsl% pe_x)
        pe_dc% ffsl% lats  = dc (i)% glate (1)
        pe_dc% ffsl% nlatx = dc (i)% nglat / 2
        IF (set_a == nproca) pe_dc% ffsl% lats  = pe_dc% glate (2)
      ELSE                            ! keep southern hemisphere on this pe
        IF(set_a == nproca) THEN
          pe_dc% ffsl% pe_x = pe_dc% mapmesh (set_b,set_a-1)
          pe_dc% ffsl% pe_n = pe_dc% mapmesh (set_b,set_a-1)
        ELSE IF(set_a == nproca-1) THEN
          pe_dc% ffsl% pe_x = pe_dc% mapmesh (set_b,set_a-1)
          pe_dc% ffsl% pe_n = pe_dc% mapmesh (set_b,set_a+1)
        ELSE
          pe_dc% ffsl% pe_x = pe_dc% mapmesh (set_b,set_a-1)
          pe_dc% ffsl% pe_n = pe_dc% mapmesh (set_b,set_a+2)
        ENDIF
        IF(set_a /= 2) THEN
          pe_dc% ffsl% pe_s = pe_dc% mapmesh (set_b,set_a-2)
        ENDIF
        pe_dc% ffsl% latn = pe_dc% glats (2)
        i = indx (pe_dc% ffsl% pe_x)
        pe_dc% ffsl% lats  = dc (i)% glate (2)
        pe_dc% ffsl% nlatx = dc (i)% nglat / 2
      ENDIF
    ENDIF
    pe_dc% ffsl% nlat = pe_dc% ffsl% lats - pe_dc% ffsl% latn + 1
    pe_dc% ffsl% lats = pe_dc% nlat + 1   - pe_dc% ffsl% lats
    pe_dc% ffsl% latn = pe_dc% nlat + 1   - pe_dc% ffsl% latn
  CONTAINS
!..............................................................................
    FUNCTION indx (pe)
    INTEGER ,INTENT(in) :: pe
    INTEGER             :: indx
      INTEGER :: i
      indx = -1
      DO i=1,SIZE(dc)
        IF(dc(i)%pe == pe) THEN
          indx = i
          EXIT
        ENDIF
      END DO
    END FUNCTION indx
  END SUBROUTINE decompose_ffsl
!------------------------------------------------------------------------------
  SUBROUTINE print_ffsl_dec (dc)
  TYPE (pe_decomposed) ,INTENT (in) :: dc (:) ! global decomposition array
    INTEGER :: i
    IF (p_pe == p_io) THEN
      IF(SIZE(dc)>1) THEN
!OSBUTLS        CALL message('', separator)
      ENDIF
      CALL message('', '')
      CALL message('','  ffsl Flux-Form Semi-Lagrangian transport scheme decomposition')
      CALL message('', '')
      CALL message('','  seta  setb    pe  pe_x  pe_n  pe_s  nlat  latn  lats nlatx')
      CALL message('', '')
      DO i=1, SIZE(dc)
        WRITE(message_text,'(10i6)')  &
          dc(i)%       set_a, &
          dc(i)%       set_b, &
          dc(i)%       pe,    &
          dc(i)% ffsl% pe_x,  &
          dc(i)% ffsl% pe_n,  &
          dc(i)% ffsl% pe_s,  &
          dc(i)% ffsl% nlat,  &
          dc(i)% ffsl% latn,  &
          dc(i)% ffsl% lats,  &
          dc(i)% ffsl% nlatx
        CALL message('', message_text)
      END DO

    ENDIF
  END SUBROUTINE print_ffsl_dec
!------------------------------------------------------------------------------
  SUBROUTINE cleanup_decomposition
    INTEGER :: i
    DO i = lbound(global_decomposition,1), ubound(global_decomposition,1)
      CALL cleanup_dc (global_decomposition(i))
    END DO
    DEALLOCATE (global_decomposition)
  END SUBROUTINE cleanup_decomposition
!------------------------------------------------------------------------------
  SUBROUTINE cleanup_dc (dc)
  TYPE (pe_decomposed) ,INTENT (inout) :: dc
    DEALLOCATE (dc% mapmesh)    
    DEALLOCATE (dc% nnp)
    DEALLOCATE (dc% nmp)
    DEALLOCATE (dc% glat)
    DEALLOCATE (dc% glon)
    IF (ASSOCIATED(dc%mask)) DEALLOCATE (dc% mask)
    DEALLOCATE (dc% lm)
    DEALLOCATE (dc% nlmp)
    DEALLOCATE (dc% nlnp)
    DEALLOCATE (dc% intr)
    DEALLOCATE (dc% np1)
    DEALLOCATE (dc% np1i)
    DEALLOCATE (dc% mymsp)
    DEALLOCATE (dc% nindex)
    DEALLOCATE (dc% sm)
    DEALLOCATE (dc% snnp)
    DEALLOCATE (dc% snn0)
  END SUBROUTINE cleanup_dc

END MODULE mo_decomposition
