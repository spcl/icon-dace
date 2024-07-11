!
! Prepares and provides cloud optical properties
!
! Remarks
!   This code makes use tabulated Mie Calculations for the RRTM band structure
!   prepared by Stefan Kinne.
!
! Origin
!   Rewrite and synthesis of ECHAM5 code, particularly rad_int.f90 
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

MODULE mo_newcld_optics

  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: finish

  USE mo_math_constants,     ONLY: pi
  USE mo_physical_constants, ONLY: rhoh2o

  USE mo_netcdf,             ONLY: nf90_nowrite, nf90_noerr
  USE mo_netcdf_parallel,    ONLY: p_nf90_open, p_nf90_close, &
    &                              p_nf90_inq_varid,          &
    &                              p_nf90_get_var

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: setup_newcld_optics, newcld_optics
  PUBLIC :: relmin, relmax, reimax, reimin


  INTEGER, PARAMETER :: &
    &  nbnds_lw   = 16,                &!< n of longwave bands
    &  nbnds_sw   = 14,                &!< n of shortwave bands
    &  n_mdl_bnds = nbnds_lw+nbnds_sw, &!< n of bands in effective radius table
    &  n_sizes    = 61                  !< n of bands in effective radius table

  REAL (wp), PARAMETER ::     &
    &  ccwmin    = 1.e-7_wp,  & !< min condensate for lw cloud opacity
    &  zkap_cont = 1.143_wp,  & !< continental (Martin et al. ) breadth param
    &  zkap_mrtm = 1.077_wp     !< maritime (Martin et al.) breadth parameter

  REAL (wp), PARAMETER, DIMENSION(16) ::rebcug = (/ &
    &  0.718_wp, 0.726_wp, 1.136_wp, 1.320_wp,      &
    &  1.505_wp, 1.290_wp, 0.911_wp, 0.949_wp,      &
    &  1.021_wp, 1.193_wp, 1.279_wp, 0.626_wp,      &
    &  0.647_wp, 0.668_wp, 0.690_wp, 0.690_wp     /)

  REAL (wp), PARAMETER, DIMENSION(16) ::rebcuh = (/ &
    &  0.0069_wp, 0.0060_wp, 0.0024_wp, 0.0004_wp,  &
    & -0.0016_wp, 0.0003_wp, 0.0043_wp, 0.0038_wp,  &
    &  0.0030_wp, 0.0013_wp, 0.0005_wp, 0.0054_wp,  &
    &  0.0052_wp, 0.0050_wp, 0.0048_wp, 0.0048_wp /)

  LOGICAL,   SAVE :: l_variable_inhoml

  REAL (wp), SAVE :: &
    &  zinpar,       & !< exponent for variable liquid-cloud inhomogeneity
    &  zinhomi         !< ice-cloud inohomogeneity factor

  REAL (wp), SAVE ::                 &
    &  wavenumber(n_mdl_bnds)      , & !< effective wavenumber for table band
    &  wavelength(n_mdl_bnds)      , & !< effective wavelength for table band
    &  re_droplet(n_sizes)         , & !< effective drop radius for table sizes
    &  re_crystal(n_sizes)         , & !< effective ice radius for table sizes
    &  reimin                      , & !< min ice effective radius (microns)
    &  reimax                      , & !< max ice effective radius (microns)
    &  relmin                      , & !< min liquid effective radius (microns)
    &  relmax                      , & !< max liquid effective radius (microns)
    &  del_rei                     , & !< ice effective radius (microns) increment
    &  del_rel                     , & !< liq effective radius (microns) increment
    &  z_ext_l(n_sizes, n_mdl_bnds), & !< tabulated extinction liquid
    &  z_ext_i(n_sizes, n_mdl_bnds), & !< tabulated extinction ice
    &  z_coa_l(n_sizes, n_mdl_bnds), & !< tabulated co-albedo liquid
    &  z_coa_i(n_sizes, n_mdl_bnds), & !< tabulated co-albedo ice
    &  z_asy_l(n_sizes, n_mdl_bnds), & !< tabulated asymmetry liquid
    &  z_asy_i(n_sizes, n_mdl_bnds)    !< tabulated asymmetry ice

CONTAINS
  !-----------------------------------------------------------------------------
  !>
  !! Sets resolution dependent parameters for cloud optics
  !
  SUBROUTINE setup_newcld_optics(data_filename)

    !> NetCDF file with RRTM Cloud Optical Properties for ECHAM6
    CHARACTER (LEN=*), INTENT(IN) :: data_filename

    INTEGER :: nf_file_id     !< id number of netcdf file
    INTEGER :: nf_var_id      !< id number of variable in netcdf file
    INTEGER :: nf_status      !< return status of netcdf function

    zinpar  = 0.10_wp
    zinhomi = 0.80_wp

    l_variable_inhoml = .FALSE.
    nf_status = p_nf90_open (TRIM(data_filename), nf90_nowrite, nf_file_id)
    !
    IF (nf_status /= nf90_noerr) THEN
      CALL finish('mo_newcld_optics/setup_newcld_optics',      &
        &         'File '//TRIM(data_filename)//' cannot be opened')
    END IF

    nf_status = p_nf90_inq_varid(nf_file_id, 'wavenumber', nf_var_id)
    nf_status = p_nf90_get_var  (nf_file_id, nf_var_id, wavenumber, (/1/),(/n_mdl_bnds/))

    nf_status = p_nf90_inq_varid(nf_file_id, 'wavelength', nf_var_id)
    nf_status = p_nf90_get_var  (nf_file_id, nf_var_id, wavelength, (/1/),(/n_mdl_bnds/))

    nf_status = p_nf90_inq_varid(nf_file_id, 're_droplet', nf_var_id)
    nf_status = p_nf90_get_var  (nf_file_id, nf_var_id, re_droplet, (/1/),(/n_sizes/))

    nf_status = p_nf90_inq_varid(nf_file_id, 're_crystal', nf_var_id)
    nf_status = p_nf90_get_var  (nf_file_id, nf_var_id, re_crystal, (/1/),(/n_sizes/))

    nf_status = p_nf90_inq_varid(nf_file_id, 'extinction_per_mass_droplet', nf_var_id)
    nf_status = p_nf90_get_var  (nf_file_id, nf_var_id, z_ext_l, (/1,1/),(/n_sizes,n_mdl_bnds/))

    nf_status = p_nf90_inq_varid(nf_file_id, 'co_albedo_droplet', nf_var_id)
    nf_status = p_nf90_get_var  (nf_file_id, nf_var_id, z_coa_l, (/1,1/),(/n_sizes,n_mdl_bnds/))

    nf_status = p_nf90_inq_varid(nf_file_id, 'asymmetry_factor_droplet', nf_var_id)
    nf_status = p_nf90_get_var  (nf_file_id, nf_var_id, z_asy_l, (/1,1/),(/n_sizes,n_mdl_bnds/))

    nf_status = p_nf90_inq_varid(nf_file_id, 'extinction_per_mass_crystal', nf_var_id)
    nf_status = p_nf90_get_var  (nf_file_id, nf_var_id, z_ext_i, (/1,1/),(/n_sizes,n_mdl_bnds/))

    nf_status = p_nf90_inq_varid(nf_file_id, 'co_albedo_crystal', nf_var_id)
    nf_status = p_nf90_get_var  (nf_file_id, nf_var_id, z_coa_i, (/1,1/),(/n_sizes,n_mdl_bnds/))

    nf_status = p_nf90_inq_varid(nf_file_id, 'asymmetry_factor_crystal', nf_var_id)
    nf_status = p_nf90_get_var  (nf_file_id, nf_var_id, z_asy_i, (/1,1/),(/n_sizes,n_mdl_bnds/))

    nf_status = p_nf90_close (nf_file_id)

    reimin = MINVAL(re_crystal)
    reimax = MAXVAL(re_crystal)
    del_rei= (re_crystal(2) - re_crystal(1))
    relmin = MINVAL(re_droplet)
    relmax = MAXVAL(re_droplet)
    del_rel= (re_droplet(2) - re_droplet(1))

  END SUBROUTINE setup_newcld_optics
  !-----------------------------------------------------------------------------
  !>
  !! Calculates cloud optical from cloud physical properties
  !!
  !
  SUBROUTINE newcld_optics(                                                 &
    & jce          ,kbdim        ,klev         ,nb_lw        ,nb_sw        ,&
    & zglac        ,zland        ,ktype        ,icldlyr      ,zt           ,&
    & zlwp         ,ziwp         ,zlwc         ,ziwc         ,zcdnc        ,&
    & zreff_liq     ,zreff_frz   ,icpl_reff    ,                            &
    & tau_lw       ,tau_sw       ,omg          ,asy                         )

    INTEGER, INTENT (IN)    ::     &
      &  jce,                      & !< actual length of grid point block
      &  kbdim,                    & !< maximum length of grid point block
      &  klev,                     & !< number of layers
      &  ktype(kbdim),             & !< type of convection
      &  nb_sw,                    & !< number of SW bands
      &  nb_lw,                    & !< number of LW bands
      &  icpl_reff,                & !< option for couplig radiation with reff
      &  icldlyr(kbdim,klev)         !< 0/1 flag for clear or cloudy layer

    REAL (wp), INTENT (IN)  ::     &
      &  zt(kbdim,klev),           & !< temperature
      &  zlwp(kbdim,klev),         & !< liquid water path
      &  ziwp(kbdim,klev),         & !< ice water path
      &  zlwc(kbdim,klev),         & !< liquid water content [g/m3]
      &  ziwc(kbdim,klev)            !< ice water content    [g/m3]

    REAL (wp), POINTER, INTENT (IN) ::     &
      &  zglac(:),                 & !< fraction of land covered by glaciers
      &  zland(:)                  !< land-sea mask. (1. = land, 0. = sea/lakes)
    REAL (wp), INTENT (IN)  ::     &
      &  zcdnc(kbdim,klev),               & !< cloud drop number concentration
      &  zreff_liq(kbdim,klev),           & !< effective radius of the liquid phase
      &  zreff_frz(kbdim,klev)              !< effective radius of the frozen phase


    REAL (wp), INTENT (OUT) ::     &
      &  tau_lw(kbdim,klev,nb_lw), & !< LW optical depth
      &  tau_sw(kbdim,nb_sw,klev), & !< SW optical depth
      &  omg(kbdim,nb_sw,klev),    & !< cloud single scattering albedo
      &  asy(kbdim,nb_sw,klev)       !< cloud asymmetry factor


    INTEGER  :: iband, jk, jl, ml1, ml2, mi1, mi2

    REAL(wp) :: ztol, ztoi, zol, zoi, zgl, zgi, wl1, wl2, wi1, wi2, zfact
    REAL(wp) :: zmsald, zmsaid

    REAL(wp) :: &
      &  zkap(kbdim),        & !< breath parameter for scaling effective radius
      &  zlwpt(kbdim),       & !< liquid water path
      &  zinhoml(kbdim),     & !< cloud inhomogeneity factor (liquid)
      &  re_droplets(kbdim), & !< effective radius of liquid distribution
      &  re_crystals(kbdim), & !< effective radius of ice distribution
      &  zscratch,           & !< effective radius of ice distribution
      &  ztau(kbdim,klev,n_mdl_bnds), & !< SW optical depth
      &  zomg(kbdim,klev,n_mdl_bnds), & !< cloud single scattering albedo
      &  zasy(kbdim,klev,n_mdl_bnds)    !< cloud asymmetry factor

    ! First check for consistency of number of LW and SW bands
    IF (nbnds_lw /= nb_lw .OR. nbnds_sw /= nb_sw) THEN
      CALL finish('mo_newcld_optics: Inconsistent number of spectral bands')
    END IF

    IF (relmin < 1.5_wp .OR. relmin > 2.5_wp ) THEN
      CALL finish('Apparently unsuccessful loading of optical tables')
    END IF
    !
    ! 1.0 Basic cloud properties
    ! --------------------------------
    IF (l_variable_inhoml) THEN
      zlwpt(1:jce) = 0.0_wp
      DO jk = 1, klev
        DO jl = 1, jce
          zlwpt(jl) = zlwpt(jl)+zlwp(jl,jk)
        END DO
      END DO
      WHERE (zlwpt(1:jce) > 1.0_wp)
        zinhoml(1:jce) = zlwpt(1:jce)**(-zinpar)
      ELSEWHERE
        zinhoml(1:jce) = 1.0_wp
      END WHERE
    ELSE
      DO jl = 1, jce
        IF(ktype(jl) .EQ. 0) THEN
          zinhoml(jl) = 0.77_wp
        ELSE
          zinhoml(jl) = 0.77_wp
        END IF
      END DO
    END IF

    IF ( icpl_reff==0 ) THEN ! Coupling of reff with radiation
      DO jl=1,jce
        zkap(jl) = zkap_cont*(zland(jl)-zglac(jl)) + zkap_mrtm*(1.0_wp-zland(jl)+zglac(jl))
      ENDDO
    END IF

    !
    ! 2.0 Cloud Optical Properties by interpolating tables in effective radius
    ! --------------------------------

    zfact = 1.0e6_wp*(3.0e-9_wp/(4.0_wp*pi*rhoh2o))**(1.0_wp/3.0_wp)
    DO jk=1,klev
      IF ( icpl_reff==1 ) THEN ! Coupling of reff with radiation
        DO jl=1,jce
            re_crystals(jl) = MAX(reimin,MIN(reimax,1.0e6_wp* zreff_frz(jl,jk)))
            re_droplets(jl) = MAX(relmin,MIN(relmax,1.0e6_wp* zreff_liq(jl,jk)))
        ENDDO
      ELSE IF ( icpl_reff==0 ) THEN
        DO jl=1,jce
          IF (icldlyr(jl,jk)==1 .AND. (zlwp(jl,jk)+ziwp(jl,jk))>ccwmin) THEN
            ! see ECHAM5 documentation (Roeckner et al, MPI report 349)            
            re_crystals(jl) = MAX(reimin ,MIN(reimax  ,83.8_wp*ziwc(jl,jk)**0.216_wp))
!mk:opt   re_crystals = MAX(20.0_wp,MIN(150.0_wp,83.8_wp*ziwc(jl,jk)**0.216_wp))
            re_droplets(jl) = MAX(relmin,MIN(relmax,zfact*zkap(jl)*(zlwc(jl,jk) &
                                       & /zcdnc(jl,jk))**(1.0_wp/3.0_wp)))
          ENDIF
        ENDDO
      ELSE
        CALL finish('Invalid option for icpl_reff in newcld_optics')
      ENDIF

!$NEC ivdep 
      DO jl=1,jce
        IF (icldlyr(jl,jk)==1 .AND. (zlwp(jl,jk)+ziwp(jl,jk))>ccwmin) THEN

          ! alternative formulation Ou and Liou (1995) as function of T as in IFS (cy38r2)
          ! re_crystals = MAX(20.0_wp, MIN(70.0_wp, 0.5_wp * &  ! limits to range of data
          !   & ( 326.3_wp                                   &
          !   & + 12.42_wp  * (zt(jl,jk)-tmelt)              &
          !   & + 0.197_wp  * (zt(jl,jk)-tmelt)**2           &
          !   & + 0.0012_wp * (zt(jl,jk)-tmelt)**3 )))

          ! optional tuning of effective crystal radius
          ! re_crystals = re_crystals * 1.25_wp

          ml1 = MAX(1,MIN(n_sizes-1,FLOOR(1.0_wp+(re_droplets(jl)-relmin)/del_rel)))
          ml2 = ml1 + 1
          wl1 = 1.0_wp - (re_droplets(jl) - (relmin + del_rel* REAL(ml1-1,wp)) )/del_rel
          wl2 = 1.0_wp - wl1

          mi1 = MAX(1,MIN(n_sizes-1,FLOOR(1.0_wp+(re_crystals(jl)-reimin)/del_rei)))
          mi2 = mi1 + 1
          wi1 = 1.0_wp - (re_crystals(jl) - (reimin + del_rei * REAL(mi1-1,wp)) )/del_rei
          wi2 = 1.0_wp - wi1

          ! Note: The following implementation is in principle intended; the first SW
          !       band overlaps with the last LW band, and the last SW band is unused
          !  ***  An open question is, however, is the usage of the coefficients
          !       in the overlapping band: currently, ztau and zomg are set for LW,
          !       whereas zasy are set for SW. ***
!$NEC unroll(nbnds_sw)
          DO iband = nbnds_lw,n_mdl_bnds-1
            ztol = zlwp(jl,jk)*(wl1*z_ext_l(ml1,iband) + wl2*z_ext_l(ml2,iband))
            ztoi = ziwp(jl,jk)*(wi1*z_ext_i(mi1,iband) + wi2*z_ext_i(mi2,iband))
            zol  = 1.0_wp - (wl1*z_coa_l(ml1,iband) + wl2*z_coa_l(ml2,iband))
            zoi  = 1.0_wp - (wi1*z_coa_i(mi1,iband) + wi2*z_coa_i(mi2,iband))
            zgl  = wl1*z_asy_l(ml1,iband) + wl2*z_asy_l(ml2,iband)
            zgi  = wi1*z_asy_i(mi1,iband) + wi2*z_asy_i(mi2,iband)

            zscratch = (ztol*zol+ztoi*zoi)
            ztau(jl,jk,iband) = ztol*zinhoml(jl) + ztoi*zinhomi
            zomg(jl,jk,iband) = zscratch/(ztol+ztoi)
            zasy(jl,jk,iband) = (ztol*zol*zgl+ztoi*zoi*zgi)/zscratch
          END DO
          !
          ! set old Cloud Optics instead of Kinne Optics for LW
          !
!$NEC unroll(nbnds_lw)
          DO iband = 1,nbnds_lw
            zmsald=0.025520637_wp+0.2854650784_wp*EXP(-0.088968393014_wp  &
                 *re_droplets(jl))
            zmsaid=(rebcuh(iband)+rebcug(iband)/re_crystals(jl))
            ztau(jl,jk,iband)  = zmsald*zlwp(jl,jk)*zinhoml(jl)          &
                 &             + zmsaid*ziwp(jl,jk)*zinhomi
            zomg(jl,jk,iband) = 0.0_wp
          END DO

        ELSE
!$NEC unroll(n_mdl_bnds-1)
          ztau(jl,jk,1:n_mdl_bnds-1)  = 0.0_wp
!$NEC unroll(n_mdl_bnds-1)
          zomg(jl,jk,1:n_mdl_bnds-1)  = 1.0_wp
!$NEC unroll(n_mdl_bnds-1)
          zasy(jl,jk,1:n_mdl_bnds-1)  = 0.0_wp
        END IF
!$NEC unroll(nbnds_lw)
        tau_lw(jl,jk,1:nbnds_lw) = ztau(jl,jk,1:nbnds_lw) * (1.0_wp - zomg(jl,jk,1:nbnds_lw))
!$NEC unroll(nbnds_sw)
        tau_sw(jl,1:nbnds_sw,jk) = ztau(jl,jk,nbnds_lw:nbnds_lw+nbnds_sw-1)
!$NEC unroll(nbnds_sw)
        omg(jl,1:nbnds_sw,jk)    = zomg(jl,jk,nbnds_lw:nbnds_lw+nbnds_sw-1)
!$NEC unroll(nbnds_sw)
        asy(jl,1:nbnds_sw,jk)    = zasy(jl,jk,nbnds_lw:nbnds_lw+nbnds_sw-1)
      END DO
    END DO

  END SUBROUTINE newcld_optics

END MODULE mo_newcld_optics
