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

! Cariolle scheme.
! documentation: cr2016_10_22_rjs

MODULE mo_lcariolle

  USE mo_kind, ONLY: wp => dp, wi => i4
  USE mo_read_interface, ONLY: read_bcast_REAL_3D, read_1D,  &
    & closeFile, openInputFile
  USE mo_physical_constants, ONLY: avo
  USE mo_fortran_tools, ONLY: set_acc_host_or_device

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_avi, t_time_interpolation, lcariolle_init_o3, &
    & l_cariolle_initialized_o3, lcariolle_init, lcariolle_do3dt 

! number of latitudes, pressure layers and months used in the Cariolle
! climatology
  INTEGER(wi), PARAMETER :: nlatx=65, nlevx=91, nmonthx=12

! The derived type t_avi is meant to create variables containing all
! variables that have to be filled in the host model and being passed to
! the Cariolle submodel.
  TYPE t_avi
    REAL(wp), POINTER, DIMENSION(:,:) :: tmprt, vmr2molm2, pres, o3_vmr
    LOGICAL, POINTER  :: lday(:)
    REAL(wp), POINTER :: cell_center_lat(:)
    LOGICAL           :: ldown
  END TYPE t_avi

! The derived type t_pvi is meant to host all variables connected to
! the Cariolle climatology of coefficients determining the ozone tendencies.
  TYPE t_pvi
    REAL(wp), DIMENSION(0:nlatx+1,nlevx,0:nmonthx+1) :: &
      & a1, a2, a3, a4, a5, a6, a7, a8
    REAL(wp) :: plev(nlevx), rlat(0:nlatx+1), delta_lat
    LOGICAL           :: l_lat_sn
  END type t_pvi

  TYPE(t_pvi), TARGET :: pvi
! The variables of type t_time_interpolation contain the weights and indices
! for a linear time interpolation of the monthly given climatology of
! coefficients of the Cariolle scheme.
  TYPE t_time_interpolation
    INTEGER(wi)          :: imonth1,imonth2
    REAL(wp)             :: weight1,weight2
  END TYPE t_time_interpolation

  LOGICAL :: l_cariolle_initialized_o3 = .FALSE.

  CHARACTER(*), PARAMETER :: modname="mo_lcariolle"

CONTAINS

!! This subroutine initializes the data structures for the
!! Cariolle interactive ozone scheme Cariolle et al.:
!! Atmos. Chem. Phys. 7, 2183 (2007), and reads all necessary coefficients.
!! documentation: cr2016_10_22_rjs
  SUBROUTINE lcariolle_init()
    INTEGER :: file_id

    CALL openInputFile(file_id, 'cariolle_coeff.nc')
    CALL read_1d(file_id, 'lat', fill_array=pvi%rlat(1:nlatx))
    CALL read_1d(file_id, 'plev', fill_array=pvi%plev)
    pvi%plev=pvi%plev*100._wp
    IF (pvi%rlat(1)<0._wp) THEN
      pvi%l_lat_sn=.TRUE.
      pvi%rlat(0)=-90._wp
      pvi%rlat(nlatx+1)=90._wp
    ELSE
      pvi%l_lat_sn=.FALSE.
      pvi%rlat(0)=90._wp
      pvi%rlat(nlatx+1)=-90._wp
    END IF
    pvi%rlat=(acos(-1._wp)/180._wp)*pvi%rlat
    pvi%delta_lat=ABS(pvi%rlat(1)-pvi%rlat(2))
    CALL read_fill_array('a1', pvi%a1)
    CALL read_fill_array('a2', pvi%a2)
    CALL read_fill_array('a3', pvi%a3)
    CALL read_fill_array('a4', pvi%a4)
    CALL read_fill_array('a5', pvi%a5)
    CALL read_fill_array('a6', pvi%a6)
    CALL read_fill_array('a7', pvi%a7)
    CALL read_fill_array('a8', pvi%a8)
    CALL closeFile(file_id)
    !$ACC ENTER DATA COPYIN(pvi, pvi%rlat, pvi%plev) &
    !$ACC   COPYIN(pvi%a1, pvi%a2, pvi%a3, pvi%a4, pvi%a5, pvi%a6, pvi%a7, pvi%a8)
  CONTAINS

    SUBROUTINE read_fill_array(aname, arr)
      CHARACTER(*), INTENT(IN) :: aname
      REAL(wp), INTENT(OUT) :: arr(:,:,0:)
      REAL(wp) :: a_3d(nmonthx,nlevx,nlatx)

      CALL read_bcast_REAL_3D(file_id, aname, fill_array=a_3d)
      arr(1:nlatx,:,1:12)=RESHAPE(a_3d,(/nlatx,nlevx,nmonthx/),ORDER=(/3,2,1/))
      arr(1:nlatx,:,0)=arr(1:nlatx,:,12)
      arr(1:nlatx,:,13)=arr(1:nlatx,:,1)
    END SUBROUTINE read_fill_array
  END SUBROUTINE lcariolle_init

!! This subroutine initializes the ozone concentration at the beginning
!! of a simulation. It has to be called outside the time loop and only once
!! for every column.
!! The subroutine can take any number of columns jcb,...,jce
!! documentation: cr2016_10_22_rjs
!!
  SUBROUTINE lcariolle_init_o3(                                   &
           & jcb,                 jce,         NCX,               &
           & nlev,                time_ip,     avi,         vmr_o3)
    INTEGER(wi),INTENT(IN) :: &
         & jcb,jce,   & !< begin, end index of column
         & NCX,       & !< first dim of fields as in calling subprogram
         & nlev         !< number of levels in column
    TYPE(t_time_interpolation), INTENT(IN) :: time_ip  !< contains linear interpolation weights
    TYPE(t_avi),INTENT(IN) :: avi                      !< derived type containing all variables
                                                       !< passed to submodel Cariolle 
    REAL(wp),INTENT(INOUT) :: vmr_o3(NCX,nlev) ! initial ozone volume mixing ratio
    INTEGER(wi)            :: ilev,ic
    REAL(wp)               :: wgt1_lat(NCX),wgt2_lat(NCX), &
                            & wgt1_p(NCX,nlev),wgt2_p(NCX,nlev)
    INTEGER(wi)            :: inmw1_lat(NCX),inmw2_lat(NCX), &
                            & iw1_p(NCX,nlev),iw2_p(NCX,nlev)
    REAL(wp)               :: wgt1,wgt2,wp1,wp2
    INTEGER(wi)            :: iw1,iw2,ip1,ip2
    REAL(wp)               :: a3_p1,a3_p2
    REAL(wp)               :: at3(0:nlatx+1,nlevx)
    
    ! calculate linear interpolation weights for latitude interpolation
    CALL lcariolle_lat_intp_li(                                      &
       & jcb,                 jce,                NCX,          &
       & avi%cell_center_lat, nlatx,              pvi%rlat,     &
       & pvi%delta_lat,       pvi%l_lat_sn,       wgt1_lat,     &
       & wgt2_lat,            inmw1_lat,          inmw2_lat     )
    ! calculate linear interpolation weights for pressure interpolation 
    CALL lcariolle_pres_intp_li(                                     &
       & jcb,                 jce,                NCX,          &
       & nlev,                avi%pres,           nlevx,        &
       & pvi%plev,            wgt1_p,             wgt2_p,       &
       & iw1_p,               iw2_p                             )
    ! interpolate coefficients with respect to time
    wp1=time_ip%weight1
    wp2=time_ip%weight2
    ip1=time_ip%imonth1
    ip2=time_ip%imonth2
    at3(:,:)=wp1*pvi%a3(:,:,ip1)+wp2*pvi%a3(:,:,ip2)
    ! latitude and pressure interpolation of at1,...,at8
    DO ilev=1,nlev
      DO ic=jcb,jce
        wgt1=wgt1_lat(ic)
        wgt2=wgt2_lat(ic)
        iw1 =inmw1_lat(ic)
        iw2 =inmw2_lat(ic)
        wp1 =wgt1_p(ic,ilev)
        wp2 =wgt2_p(ic,ilev)
        ip1 =iw1_p(ic,ilev)
        ip2 =iw2_p(ic,ilev)
    ! latitude interpolation for pressure level 1 (iw1_p)
        a3_p1 = wgt1*at3(iw1,ip1) + &
              & wgt2*at3(iw2,ip1)
    ! latitude interpolation for pressure level 2 (iw2_p)
        a3_p2 = wgt1*at3(iw1,ip2) + &
              & wgt2*at3(iw2,ip2)
    ! pressure level interpolation
        vmr_o3(ic,ilev)=wp1*a3_p1+wp2*a3_p2
      END DO
    END DO
  END SUBROUTINE lcariolle_init_o3

!! This subroutine calculates interpolation weights for a linear
!! interpolation with respect to pressure layers.
!! documentation: cr2016_10_22_rjs
!!
  SUBROUTINE lcariolle_pres_intp_li (                                 &
             & jcb,              jce,                 NCX,            &
             & nlev,             plev,                nlev_clim,      &
             & plev_clim,        wgt1_p,              wgt2_p,         &
             & inmw1_p,          inmw2_p,                             &
             & lacc)
    INTEGER(wi),INTENT(IN)        :: &
         & jcb,jce,  & !< begin, end index of column
         & NCX,      & !< first dim of fields as in calling subprogram
         & nlev        !< number of levels in column
    REAL(wp),INTENT(IN)           :: plev(NCX,nlev) !< pressure to interpolate on
    INTEGER(wi),INTENT(IN)        :: nlev_clim      !< number of pressures in climat.
    REAL(wp),INTENT(IN)           :: plev_clim(nlev_clim) !< climatological pressures
    REAL(wp),INTENT(OUT)          :: wgt1_p(NCX,nlev),wgt2_p(NCX,nlev) !< intp. weights
    INTEGER(wi),INTENT(OUT)       :: inmw1_p(NCX,nlev),inmw2_p(NCX,nlev) !< indices
    INTEGER(wi)                   :: ic,ilev,ii,ilev_low
    LOGICAL, OPTIONAL, INTENT(in) :: lacc !< GPU flag
    LOGICAL                       :: lzacc ! non-optional version of lacc

    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR SEQ ASYNC(1) PRIVATE(ilev_low, ii) IF(lzacc)
    DO ic=jcb,jce
      ilev_low=1
      DO ilev=1,nlev
        ii=ilev_low
        DO
          ! Pressures are assumed to be in irregular distances in the climatology
          IF (ii>nlev_clim) THEN
            inmw1_p(ic,ilev)=nlev_clim
            inmw2_p(ic,ilev)=nlev_clim
            wgt1_p(ic,ilev)=0.5_wp
            wgt2_p(ic,ilev)=0.5_wp
            EXIT
          END IF
          IF(plev(ic,ilev).GE.plev_clim(ii)) THEN
            ii=ii+1
          ELSE
            IF(ii==1) THEN
              inmw1_p(ic,ilev)=1
              inmw2_p(ic,ilev)=1
              wgt1_p(ic,ilev)=0.5_wp
              wgt2_p(ic,ilev)=0.5_wp
            ELSE
              inmw1_p(ic,ilev)=ii-1
              inmw2_p(ic,ilev)=ii
              wgt1_p(ic,ilev)=(plev_clim(ii)-plev(ic,ilev))/ &
                             &(plev_clim(ii)-plev_clim(ii-1))
              wgt2_p(ic,ilev)=1._wp-wgt1_p(ic,ilev)
            END IF
            EXIT
          END IF
        END DO
        ilev_low=ii
      END DO
    END DO
    !$ACC END PARALLEL LOOP
  END SUBROUTINE lcariolle_pres_intp_li

!! Interpolation weights for linear interpolation from a
!! regular grid between N pole and S pole. The N pole and S pole have to
!! be included, but the distance to the first (last) point of the regular
!! grid can be different. 
!! documentation: cr2016_10_22_rjs
!!
  SUBROUTINE lcariolle_lat_intp_li(                                         &
             & jcb,                jce,                   NCX,              &
             & cell_center_lat,    n_lat_clim,            r_lat_clim,       &
             & r_delta_lat_clim,   l_lat_clim_sn,         wgt1_lat,         &
             & wgt2_lat,           inmw1_lat,             inmw2_lat,        &
             & lacc)
    INTEGER(wi),INTENT(IN)        :: jcb, jce, & !< first, last point to interpolate on
                                   & NCX         !< length of array in calling subprog.
    REAL(wp),INTENT(IN)           :: cell_center_lat(NCX) !< latitudes in radiant
    INTEGER(wi),INTENT(IN)        :: n_lat_clim           !< number of lats in climatol.
    REAL(wp),INTENT(IN)           :: r_lat_clim(0:n_lat_clim+1), & !< lats in climatology
                                                               !< in radiant
                                   & r_delta_lat_clim     !< spacing in climatology
    LOGICAL,INTENT(IN)            :: l_lat_clim_sn        !< .true. if order is S->N,
                                                          !< .false. otherwise
    REAL(wp),INTENT(OUT)          :: wgt1_lat(NCX),wgt2_lat(NCX)   !< interpolation weights
    INTEGER(wi),INTENT(OUT)       :: inmw1_lat(NCX),inmw2_lat(NCX) !< corresponding indices 
    INTEGER(wi)                   :: n_order
    REAL(wp)                      :: r_delta_lat_clim_i
    INTEGER                       :: ic
    LOGICAL, OPTIONAL, INTENT(in) :: lacc !< GPU flag
    LOGICAL                       :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)
    
    n_order = MERGE(1, -1, l_lat_clim_sn)
    r_delta_lat_clim_i=1._wp/r_delta_lat_clim

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1) IF(lzacc)
    DO ic=jcb,jce
       inmw1_lat(ic)=MAX(INT(n_order*(cell_center_lat(ic)-r_lat_clim(1))*r_delta_lat_clim_i+1._wp),0)
       inmw2_lat(ic)=inmw1_lat(ic)+1
       wgt2_lat(ic)=n_order*(cell_center_lat(ic)-r_lat_clim(inmw1_lat(ic)))*r_delta_lat_clim_i
       wgt1_lat(ic)=1.0_wp-wgt2_lat(ic)
    END DO
    !$ACC END PARALLEL LOOP
  END SUBROUTINE lcariolle_lat_intp_li

!! This subroutine calculates the ozone tendency according to equation (1)
!! of Cariolle et al.: Atmos. Chem. Phys. 7, 2183 (2007).
!! The subroutine can take any number of columns jcb,...,jce, but the columns must
!! comprise the full overhead ozone column.
!! documentation: cr2016_10_22_rjs
!!
  SUBROUTINE lcariolle_do3dt(                                    &
             & jcb,             jce,          NCX,               &
             & nlev,            time_ip,      avi,          do3dt)
    INTEGER(wi),INTENT(IN)                 :: &
                                            & jcb,jce, & !< begin, end index of column
                                            & NCX,     & !< first dim of fields as in calling subprogram
                                            & nlev       !< number of levels in column
    TYPE(t_time_interpolation), INTENT(IN) :: time_ip    !< contains linear interpolation weights

    TYPE(t_avi),INTENT(IN)    :: avi                     !< derived type containing all variables
                                                         !< passed to submodel Cariolle
    REAL(wp),   INTENT(INOUT) :: do3dt(NCX,nlev)         !< tendency of ozone VMR per second
    
    INTEGER(wi)            :: ilev,ic
    REAL(wp)               :: o3_column(NCX,nlev)        !< overhead ozone column 
    REAL(wp)               :: wgt1_lat(NCX),wgt2_lat(NCX), &
                            & wgt1_p(NCX,nlev),wgt2_p(NCX,nlev)
    INTEGER(wi)            :: inmw1_lat(NCX),inmw2_lat(NCX), &
                            & iw1_p(NCX,nlev),iw2_p(NCX,nlev)
    REAL(wp)               :: wgt1,wgt2,wp1,wp2
    INTEGER(wi)            :: iw1,iw2,ip1,ip2
    REAL(wp)               :: a1_p1,a1_p2,a2_p1,a2_p2,a3_p1,a3_p2,a4_p1,a4_p2, &
                            & a5_p1,a5_p2,a6_p1,a6_p2,a7_p1,a7_p2,a8_p1,a8_p2
    REAL(wp)               :: al1(NCX,nlev), al2(NCX,nlev), &
                            & al3(NCX,nlev), al4(NCX,nlev), &
                            & al5(NCX,nlev), al6(NCX,nlev), &
                            & al7(NCX,nlev), al8(NCX,nlev)
    REAL(wp)               :: at1(0:nlatx+1,nlevx), at2(0:nlatx+1,nlevx), &
                            & at3(0:nlatx+1,nlevx), at4(0:nlatx+1,nlevx), &
                            & at5(0:nlatx+1,nlevx), at6(0:nlatx+1,nlevx), &
                            & at7(0:nlatx+1,nlevx), at8(0:nlatx+1,nlevx)
    INTEGER                :: ilevx, ilatx

    !$ACC DATA PRESENT(do3dt, avi, avi%o3_vmr, avi%vmr2molm2, avi%ldown, avi%cell_center_lat, avi%pres, avi%tmprt) &
    !$ACC   CREATE(o3_column) &
    !$ACC   CREATE(wgt1_lat, wgt2_lat, inmw1_lat, inmw2_lat) &
    !$ACC   CREATE(at1, at2, at3, at4, at5, at6, at7, at8) &
    !$ACC   CREATE(wgt1_p, wgt2_p, iw1_p, iw2_p) &
    !$ACC   CREATE(al1, al2, al3, al4, al5, al6, al7, al8)

    ! calculate overhead ozone column for each model layer
    CALL lcariolle_o3_column(       &
       & jcb, jce, NCX, nlev,       &
       & o3_vmr    = avi%o3_vmr,    &
       & vmr2molm2 = avi%vmr2molm2, &
       & ldown     = avi%ldown,     &
       & o3_column = o3_column      )
    !
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
    DO ilev=1,nlev
      DO ic=jcb,jce
        o3_column(ic,ilev)=o3_column(ic,ilev)*avo*1.e-4_wp
      END DO
    END DO
    !$ACC END PARALLEL LOOP
   
    ! calculate linear interpolation weights for latitude interpolation
    CALL lcariolle_lat_intp_li(                  &
       & jcb, jce, NCX,                          &
       & cell_center_lat  = avi%cell_center_lat, &
       & n_lat_clim       = nlatx,               &
       & r_lat_clim       = pvi%rlat,            &
       & r_delta_lat_clim = pvi%delta_lat,       &
       & l_lat_clim_sn    = pvi%l_lat_sn,        &
       & wgt1_lat         = wgt1_lat,            &
       & wgt2_lat         = wgt2_lat,            &
       & inmw1_lat        = inmw1_lat,           &
       & inmw2_lat        = inmw2_lat,           &
       & lacc             = .TRUE.               )
    
    ! calculate linear interpolation weights for pressure interpolation
    CALL lcariolle_pres_intp_li( &
       & jcb, jce, NCX, nlev,    &
       & plev      = avi%pres,   &
       & nlev_clim = nlevx,      &
       & plev_clim = pvi%plev,   &
       & wgt1_p    = wgt1_p,     &
       & wgt2_p    = wgt2_p,     &
       & inmw1_p   = iw1_p,      &
       & inmw2_p   = iw2_p,      &
       & lacc      = .TRUE.      )

    ! interpolate coefficients with respect to time
    wp1=time_ip%weight1
    wp2=time_ip%weight2
    ip1=time_ip%imonth1
    ip2=time_ip%imonth2
    
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
    DO ilevx=1,nlevx
      DO ilatx=0,nlatx+1
        at1(ilatx,ilevx)=wp1*pvi%a1(ilatx,ilevx,ip1)+wp2*pvi%a1(ilatx,ilevx,ip2)
        at2(ilatx,ilevx)=wp1*pvi%a2(ilatx,ilevx,ip1)+wp2*pvi%a2(ilatx,ilevx,ip2)
        at3(ilatx,ilevx)=wp1*pvi%a3(ilatx,ilevx,ip1)+wp2*pvi%a3(ilatx,ilevx,ip2)
        at4(ilatx,ilevx)=wp1*pvi%a4(ilatx,ilevx,ip1)+wp2*pvi%a4(ilatx,ilevx,ip2)
        at5(ilatx,ilevx)=wp1*pvi%a5(ilatx,ilevx,ip1)+wp2*pvi%a5(ilatx,ilevx,ip2)
        at6(ilatx,ilevx)=wp1*pvi%a6(ilatx,ilevx,ip1)+wp2*pvi%a6(ilatx,ilevx,ip2)
        at7(ilatx,ilevx)=wp1*pvi%a7(ilatx,ilevx,ip1)+wp2*pvi%a7(ilatx,ilevx,ip2)
        at8(ilatx,ilevx)=wp1*pvi%a8(ilatx,ilevx,ip1)+wp2*pvi%a8(ilatx,ilevx,ip2)
      END DO
    END DO 
    !$ACC END PARALLEL LOOP
    
    ! latitude and pressure interpolation of at1,...,at8
    !$ACC PARALLEL DEFAULT(PRESENT) PRIVATE(wgt1, wgt2, iw1, iw2, wp1, wp2, ip1, ip2) &
    !$ACC   PRIVATE(a1_p1, a2_p1, a3_p1, a4_p1, a5_p1, a6_p1, a7_p1, a8_p1) &
    !$ACC   PRIVATE(a1_p2, a2_p2, a3_p2, a4_p2, a5_p2, a6_p2, a7_p2, a8_p2) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO ilev=1,nlev
      DO ic=jcb,jce
        wgt1=wgt1_lat(ic)
        wgt2=wgt2_lat(ic)
        iw1 =inmw1_lat(ic)
        iw2 =inmw2_lat(ic)
        wp1 =wgt1_p(ic,ilev)
        wp2 =wgt2_p(ic,ilev)
        ip1 =iw1_p(ic,ilev)
        ip2 =iw2_p(ic,ilev)
        ! latitude interpolation for pressure level 1 (iw1_p)
        a1_p1 = wgt1*at1(iw1,ip1) + &
              & wgt2*at1(iw2,ip1)
        a2_p1 = wgt1*at2(iw1,ip1) + &
              & wgt2*at2(iw2,ip1)
        a3_p1 = wgt1*at3(iw1,ip1) + &
              & wgt2*at3(iw2,ip1)
        a4_p1 = wgt1*at4(iw1,ip1) + &
              & wgt2*at4(iw2,ip1)
        a5_p1 = wgt1*at5(iw1,ip1) + &
              & wgt2*at5(iw2,ip1)
        a6_p1 = wgt1*at6(iw1,ip1) + &
              & wgt2*at6(iw2,ip1)
        a7_p1 = wgt1*at7(iw1,ip1) + &
              & wgt2*at7(iw2,ip1)
        a8_p1 = wgt1*at8(iw1,ip1) + &
             & wgt2*at8(iw2,ip1)
    ! latitude interpolation for pressure level 2 (iw2_p)
        a1_p2 = wgt1*at1(iw1,ip2) + &
              & wgt2*at1(iw2,ip2)
        a2_p2 = wgt1*at2(iw1,ip2) + &
              & wgt2*at2(iw2,ip2)
        a3_p2 = wgt1*at3(iw1,ip2) + &
              & wgt2*at3(iw2,ip2)
        a4_p2 = wgt1*at4(iw1,ip2) + &
              & wgt2*at4(iw2,ip2)
        a5_p2 = wgt1*at5(iw1,ip2) + &
              & wgt2*at5(iw2,ip2)
        a6_p2 = wgt1*at6(iw1,ip2) + &
              & wgt2*at6(iw2,ip2)
        a7_p2 = wgt1*at7(iw1,ip2) + &
              & wgt2*at7(iw2,ip2)
        a8_p2 = wgt1*at8(iw1,ip2) + &
              & wgt2*at8(iw2,ip2)
    ! pressure level interpolation
        al1(ic,ilev) = wp1*a1_p1 + wp2*a1_p2
        al2(ic,ilev) = wp1*a2_p1 + wp2*a2_p2
        al3(ic,ilev) = wp1*a3_p1 + wp2*a3_p2
        al4(ic,ilev) = wp1*a4_p1 + wp2*a4_p2
        al5(ic,ilev) = wp1*a5_p1 + wp2*a5_p2
        al6(ic,ilev) = wp1*a6_p1 + wp2*a6_p2
        al7(ic,ilev) = wp1*a7_p1 + wp2*a7_p2
        al8(ic,ilev) = wp1*a8_p1 + wp2*a8_p2
      END DO
    END DO
    !$ACC END PARALLEL
    
    ! calculate equation (1) of Cariolle et al., Atmos. Chem. Phys. 7, 2183 (2007).
    ! first step: all terms except the A_8 term for polar stratospheric clouds
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
    DO ilev=1,nlev
      DO ic=jcb,jce
        do3dt(ic,ilev) = al1(ic,ilev) + &
                       & al2(ic,ilev)*(avi%o3_vmr(ic,ilev)-al3(ic,ilev)) + &
                       & al4(ic,ilev)*(avi%tmprt(ic,ilev)-al5(ic,ilev)) + &
                       & al6(ic,ilev)*(o3_column(ic,ilev)-al7(ic,ilev))
      END DO
    END DO
    !$ACC END PARALLEL LOOP
    
    ! Add A_8 in case of polar stratospheric clouds and daylight for
    ! additional ozone destruction
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
    DO ilev=1,nlev
      DO ic=jcb,jce
        IF (avi%lday(ic).AND.avi%tmprt(ic,ilev)<=195._wp) THEN
          do3dt(ic,ilev) = do3dt(ic,ilev) + al8(ic,ilev)*avi%o3_vmr(ic,ilev)
        END IF
      END DO
    END DO
    !$ACC END PARALLEL LOOP
    !$ACC WAIT
    !$ACC END DATA
  END SUBROUTINE lcariolle_do3dt

!! This subroutine calculates the overhead ozone column in mole/m^2
!! The subroutine can take any number of columns jcb,...,jce, but the columns must
!! comprise the full ozone column.
!! documentation: cr2016_10_22_rjs
!!
  SUBROUTINE lcariolle_o3_column(                    &
           & jcb,          jce,           NCX,       &
           & nlev,         o3_vmr,        vmr2molm2, &
           & ldown,        o3_column                 )
    INTEGER(wi),INTENT(IN)  :: &
         & jcb,jce,  & !< begin, end index of column
         & NCX,      & !< first dim of fields as in calling subprogram
         & nlev        !< number of levels in column
    REAL(wp),INTENT(IN)     :: o3_vmr(NCX,nlev)    !< ozone VMR 
    REAL(wp),INTENT(IN)     :: vmr2molm2(NCX,nlev) !< conversion factor from VMR to mole/m^2
    LOGICAL,INTENT(IN)      :: ldown               !< .true. if layers are counted from top
                                                   !< to bottom, .false. otherwise
    REAL(wp),INTENT(INOUT)  :: o3_column(NCX,nlev) !< overhead ozone column (mole/m^2)
    INTEGER                 :: ic
    INTEGER(wi)             :: ii,ilev,ilev0,incr
    
    ilev = MERGE(1, nlev, ldown)
    incr = MERGE(1, -1, ldown)
    !
    !$ACC DATA PRESENT(o3_vmr, vmr2molm2, o3_column)
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO ic = jcb,jce
      o3_column(ic,ilev) = 0.5_wp * o3_vmr(ic,ilev) * vmr2molm2(ic,ilev)
    END DO
    !$ACC END PARALLEL LOOP

    ilev0 = ilev
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP SEQ
    DO ii=1,nlev-1
      !$ACC LOOP GANG VECTOR PRIVATE(ilev)
      DO ic = jcb,jce
        ilev = ilev0 + ii*incr  
        o3_column(ic,ilev) = o3_column(ic,ilev-incr) +                   &
              & 0.5_wp * o3_vmr(ic,ilev-incr) * vmr2molm2(ic,ilev-incr) +&
              & 0.5_wp * o3_vmr(ic,ilev) * vmr2molm2(ic,ilev)
      END DO
    END DO
    !$ACC END PARALLEL
    !$ACC END DATA
  END SUBROUTINE lcariolle_o3_column

END MODULE mo_lcariolle
