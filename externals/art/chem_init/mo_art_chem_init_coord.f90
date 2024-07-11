!
! mo_chem_init_coord
! This module is based on mo_vertical_coord_table
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

MODULE mo_art_chem_init_coord

  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: finish
  USE mo_impl_constants,     ONLY: success, max_char_length
  USE mo_physical_constants, ONLY: grav, rcpd, rd, p0sl_bg
! ART
  USE mo_art_chem_init_types,ONLY: t_art_chem_init_coord

  IMPLICIT NONE

  PRIVATE

  ! SUBROUTINES
  PUBLIC :: alloc_vct_chem_init 
  PUBLIC :: init_vct_chem_init
  PUBLIC :: half_level_pressure_chem_init
  PUBLIC :: auxhyb_chem_init
  PUBLIC :: full_level_pressure_chem_init
  PUBLIC :: geopot_chem_init
CONTAINS


  SUBROUTINE alloc_vct_chem_init(coord,nlev)
    !<
    ! SUBROUTINE alloc_vct_chem_init                   
    ! Part of Module: mo_chem_init_coord
    ! based on mo_vertical_coord_table
    !
    ! Author: Jennifer Schroeter, KIT
    ! Initial Release: 2015-11-15                
    ! Modifications:
    !>

    INTEGER, INTENT(IN) :: nlev
    TYPE(t_art_chem_init_coord), INTENT (inout) ::  &
      &  coord
    INTEGER :: ist
    INTEGER :: nlevp1

    CHARACTER(len=max_char_length),PARAMETER :: routine = &
         & 'mo_vertical_coord_table:alloc_vct'

    !-----------------------------------------------------------------------
    !BOC

    nlevp1 = nlev+1

    ALLOCATE (coord%ralpha_chem_init(nlev), STAT=ist)
    IF(ist/=success)THEN
      CALL finish (TRIM(routine), ' allocation of coord%ralpha failed')
    ENDIF

    ALLOCATE (coord%rlnpr_chem_init(nlev), STAT=ist)
    IF(ist/=success)THEN
      CALL finish (TRIM(routine), ' allocation of rlnpr failed')
    ENDIF

    ALLOCATE (coord%delpr_chem_init(nlev), STAT=ist)
    IF(ist/=success)THEN
      CALL finish (TRIM(routine), ' allocation of rdelpr failed')
    ENDIF

    ALLOCATE (coord%rdelpr_chem_init(nlev), STAT=ist)
    IF(ist/=success)THEN
      CALL finish (TRIM(routine), ' allocation of rdelpr failed')
    ENDIF

  END SUBROUTINE alloc_vct_chem_init


  !>
  !!  Initializes constants for vertical coordinate calculations.
  !!
  !!  Method:
  !!    Compute loop indices and surface-pressure independent
  !!    variables associated with the vertical finite-difference scheme.
  !!    Output is in module *mo_hyb*
  !!
  !! @par Revision History
  !!  A. J. Simmons, ECMWF, November 1981, original source
  !!  L. Kornblueh, MPI, May 1998, f90 rewrite
  !!  U. Schulzweida, MPI, May 1998, f90 rewrite
  !!  A. Rhodin, MPI, Jan 1999, subroutine inihyb -> module mo_hyb
  !!  H. Wan, MPI-M, Feb 2006, rewrite: "goto" removed
  !!  H. Wan, MPI-M, Aug 2007, new name: init_hyb_params
  !!  A. Gassmann, MPI-M, (2008-04-23), change coord%apzero to 10^5Pa
  !! @par
  !!  for more details see file AUTHORS
  !!



    SUBROUTINE init_vct_chem_init(coord,nlev)

    !<
    ! SUBROUTINE init_vct_chem_init                   
    ! Part of Module: mo_chem_init_coord
    ! based on mo_vertical_coord_table
    !
    ! Author: Jennifer Schroeter, KIT
    ! Initial Release: 2015-11-15                
    ! Modifications:
    !>


    INTEGER, INTENT(IN) :: nlev
    TYPE(t_art_chem_init_coord), INTENT(inout) :: &
      &  coord

    !  Local scalars:
    REAL(wp) :: za, zb, zetam, zetap, zp, zp0icao, zpp, zrd, zs, zsm
    INTEGER  :: ilev, ilevp1, iplev, iplvp1, is, ism, ist, &
      &         jk, jlev, nvclev

    !  Intrinsic functions
    INTRINSIC EXP, LOG

    !-----------------------------------------------------------------------
    !BOC

    !  Executable statements

    !-- 1. Initialize variables
    nvclev = nlev+1
!ag    coord%apzero    = 101325._wp ! changed for NCAR summer colloquium!
!sv    coord%apzero    = 100000._wp ! this was activated before adding the following line
    coord%apzero    = p0sl_bg
    zrd       = rd
    coord%ralpha_chem_init(1) = zrd*LOG(2._wp)
    coord%rlnpr_chem_init(1)  = 2._wp*coord%ralpha_chem_init(1)
    ilev      = nlev
    ilevp1    = ilev + 1
    coord%nlevm1    = ilev - 1
    iplev     = 0
    iplvp1    = 1
    is        = nvclev + ilevp1
    ism       = is - 1
    zpp       = coord%vct_chem_init(1)
    zsm       = coord%vct_chem_init(is)

    coord%t0icao  = 288._wp
    coord%tsticao = 216.5_wp
    zp0icao = p0sl_bg
    coord%rdlnp0i = rd*LOG(zp0icao)
    coord%rdtstic = rd*coord%tsticao
    coord%alrrdic = 0.0065_wp/grav
    coord%rdt0ral = coord%t0icao/coord%alrrdic
    coord%rdlnpti = coord%rdlnp0i + (LOG(coord%tsticao/coord%t0icao))/coord%alrrdic
    coord%ptricao = EXP(coord%rdlnpti/rd)
    coord%gsticao = coord%tsticao*(coord%rdlnpti-1._wp/coord%alrrdic)

    zb      = coord%vct_chem_init(nvclev+iplvp1+1)




    !-- 2. Calculate pressure-level values

    DO WHILE ( zb == 0._wp )

      iplev  = iplvp1
      iplvp1 = iplev + 1
      IF (iplvp1==ilevp1) EXIT    ! if all levels are pressure levels

      zp            = zpp
      zpp           = coord%vct_chem_init(iplvp1)
      coord%delpr_chem_init(iplev)  = zpp - zp
      coord%rdelpr_chem_init(iplev) = 1._wp/coord%delpr_chem_init(iplev)

      IF ( iplev>1 ) THEN
        coord%rlnpr_chem_init(iplev)  = zrd*LOG(zpp/zp)
        coord%ralpha_chem_init(iplev) = zrd - zp*coord%rlnpr_chem_init(iplev)  &
                       &                / coord%delpr_chem_init(iplev)
      END IF

      zb            = coord%vct_chem_init(nvclev+iplvp1+1)

    ENDDO


    IF (iplvp1/=ilevp1) THEN   ! All levels are not pressure-levels

      coord%nplev  = iplev
      coord%nplvp1 = iplvp1
      coord%nplvp2 = iplvp1 + 1
      IF (iplev==0) THEN
        coord%nplvpa = 2
      ELSE
        coord%nplvpa = iplvp1
      END IF

      !-- 3. Calculate sigma-level values
      !coord%nplev  = iplev
      !coord%nplvp1 = iplvp1
      !coord%nplvp2 = iplvp1 + 1
      za = coord%vct_chem_init(ism-nvclev)

      DO WHILE ( za == 0._wp )

        is  = ism
        ism = is - 1
        ist = is - nvclev
        zs  = zsm
        zsm = coord%vct_chem_init(is)
        IF (ist==1) THEN
          coord%nlmsgl = 0
          coord%nlmslp = 1
          coord%nlmsla = 2
          EXIT
        ELSE
          coord%rlnpr_chem_init(ist)  = zrd*LOG(zs/zsm)
          coord%ralpha_chem_init(ist) = zrd - zsm*coord%rlnpr_chem_init(ist)/(zs-zsm)
        END IF
        za = coord%vct_chem_init(ism-nvclev)


      END DO

      IF (za>0._wp) THEN
        coord%nlmsgl = ism - nvclev
        coord%nlmslp = coord%nlmsgl + 1
        coord%nlmsla = coord%nlmslp
      END IF

    ENDIF  ! If all levels are not pressure-levels

    !-- 5. Compute full level values of the hybrid coordinate

    zetam    = (coord%vct_chem_init(1))/coord%apzero + (coord%vct_chem_init(nvclev+1))

    DO jlev = 1, nlev
      zetap         = coord%vct_chem_init(jlev+1)/coord%apzero   &
          &            + (coord%vct_chem_init(nvclev+1+jlev))

      zetam = zetap
    END DO

  END SUBROUTINE init_vct_chem_init

  !>
  !! Calculate half-level pressures at all model levels
  !! for a given surface pressure.
  !!
  !! @par Method
  !!  Calculations are performed separately for pressure,
  !!  hybrid and sigma levels.
  !!
  !! @par Arguments
  !!   *ps*        surface pressure.
  !!   *kdimp*     first dimension of 2-d array *ph.*
  !!   *klen*      number of points for which calculation is
  !!               performed.
  !!   *ph*        computed half-level pressures.
  !!
  !! @par Parameters
  !!  Required constants are obtained from module *mo_hyb*.
  !!  The latter must have been initialized by a call of
  !!  subroutine init_vertical_coord.
  !!
  !! @par Results
  !!  Results are computed for *klen* consecutive points at
  !!  each model half level.
  !!
  !! @see
  !!  External documentation of the model equations and the
  !!  organization of the vertical calculation.
  !!
  !! @par Revision History
  !!    A. J. Simmons, ECMWF, November 1981, original source
  !!    L. Kornblueh, MPI, May 1998, f90 rewrite
  !!    U. Schulzweida, MPI, May 1998, f90 rewrite
  !!    H. Wan, MPI, Feb 2006, adapted for ICOHDC
  !!    H. Wan, MPI, Jan 2010, renamed the interface
  !!

  SUBROUTINE half_level_pressure_chem_init(coord, ps,kdimp,klen,nlev_in, ph)
    !<
    ! SUBROUTINE half_level_pressure_chem_init                   
    ! Part of Module: mo_chem_init_coord
    ! based on mo_vertical_coord_table
    !
    ! Author: Jennifer Schroeter, KIT
    ! Initial Release: 2015-11-15                
    ! Modifications:
    !>
    TYPE(t_art_chem_init_coord), INTENT(in) ::  &
      &  coord
    INTEGER ,INTENT(in)  :: kdimp
    REAL(wp),INTENT(in)  :: ps(kdimp)   !< surface pressure
    INTEGER ,INTENT(in)  :: klen, nlev_in

    REAL(wp),INTENT(inout) :: ph(kdimp,nlev_in+1) !< half-level pressure

    REAL(wp):: zb, zp
    INTEGER :: jk, jl, nvclev, nlevp1

    nvclev = nlev_in+1
    nlevp1 = nvclev

!    DO jk = 1, coord%nplvp1
!     DO jl = 1, klen
!        ph(jl,jk) = 100000*coord%vct_chem_init(jk) +ps(jl)* coord%vct_chem_init(jk+nvclev)
!        ENDDO
!        ENDDO



    ! Transfer pressure level values
        
    DO jk = 1, coord%nplvp1
      zp = coord%vct_chem_init(jk)

      DO jl = 1, klen
        ph(jl,jk) = zp
        
      END DO
    END DO

    ! Compute hybrid level values

    DO jk = coord%nplvp2, coord%nlmsgl
      zp = coord%vct_chem_init(jk)
      zb = coord%vct_chem_init(jk+nvclev)
      DO jl = 1, klen
        ph(jl,jk) = zp + zb*ps(jl)
      END DO
!      print*, "index:", jk, "hlp:", ph(1,jk)
    END DO

    ! Compute sigma-level values

    DO jk = coord%nlmslp, nlevp1
      zb = coord%vct_chem_init(jk+nvclev)
      DO jl = 1, klen
        ph(jl,jk) = zb*ps(jl)
      END DO
    END DO

  END SUBROUTINE half_level_pressure_chem_init


  !>
  !! Calculate full-level pressures for all vertical layers
  !! Method: Simmons and Burridge (Mon.Wea.Rev.,1981,p761,Eqn.(3.18))
  !!
  !! @par Parameters
  !!    *pres_i*    half-level pressure values.
  !!    *pres_m*    computed full-level pressure values.
  !!    *kdimp*     first dimension of 2-d array *pres_i*
  !!    *klen*      number of points for which calculation is
  !!                performed.
  !!
  !!  Required constants are obtained from module *mo_hyb*.
  !!  The latter must have been initialiazed
  !!  by a call of subroutine *inihyb*.
  !!
  !! @par Revision History
  !!    H. Wan, MPI, 2006-08-17
  !!
  SUBROUTINE full_level_pressure_chem_init(coord, pres_i, kdimp, klen, nlev_in, pres_m)
    !<
    ! SUBROUTINE full_level_pressure_chem_init                   
    ! Part of Module: mo_chem_init_coord
    ! based on mo_vertical_coord_table
    !
    ! Author: Jennifer Schroeter, KIT
    ! Initial Release: 2015-11-15                
    ! Modifications:
    !>
    TYPE(t_art_chem_init_coord), INTENT(in) :: &
      &  coord
    INTEGER ,INTENT(in) :: kdimp, klen, nlev_in    !< dimension parameters
    REAL(wp),INTENT(in) :: pres_i(kdimp,nlev_in+1) !< half-level pressure

    REAL(wp),INTENT(inout) :: pres_m(kdimp,nlev_in) !< full(/mid)-level pressure

    REAL(wp):: ztmp, zpres_i_top_min
    INTEGER :: jk, jl, ikp, ik_top, nlev

    !-----

    nlev = nlev_in
    zpres_i_top_min = coord%vct_chem_init(1)
    IF ( zpres_i_top_min > 0._wp ) THEN
      ik_top = 1
    ELSE
      ik_top = 2
      pres_m(1:klen,1) = pres_i(1:klen,2)*0.5_wp
    END IF

    DO jk = ik_top, nlev
       ikp = jk+1
       DO jl = 1, klen
         ztmp = ( pres_i(jl,ikp)*LOG(pres_i(jl,ikp))   &
         &       -pres_i(jl,jk )*LOG(pres_i(jl,jk )) ) &
         &     /( pres_i(jl,ikp)-pres_i(jl,jk) )
         pres_m(jl,jk) = EXP(ztmp-1._wp)
       END DO
    END DO

  END SUBROUTINE full_level_pressure_chem_init
  !-------------------------------------------------------------------------
  !>
  !! Calculates auxiliary variables connected with
  !! the vertical finite-difference scheme.
  !!
  !! @par Arguments
  !!
  !!   *ph*          *specified half-level pressures.
  !!   *kdim*        *first dimension of 2-d arrays *pdelp,*
  !!                  *plnpr,* *palpha,* and *ph.*
  !!   *klen*        *number of points for which calculation is performed.
  !!   *pdelp*       *computed pressure difference across layers.
  !!   *prdelp*      *reciprocal of *pdelp.*
  !!   *plnpr*       *computed logarithm of ratio of pressures.
  !!   *palpha*      *computed alphas for integration of the
  !!                  hydrostatic equation and related terms.
  !!
  !!  Required constants are obtained from modules
  !!  *mo_constants* and *mo_hyb*. The latter should have been
  !!  initialized by a call of subroutine *inihyb*.
  !!
  !!  Results are computed for *klen* consecutive points for
  !!  all required levels.
  !!
  !!  Calculations are performed separately for pressure,
  !!  hybrid and sigma levels.
  !!
  !!  External documentation of the model equations and the
  !!  organization of the vertical calculation can be found
  !!  in MPI-M technical report No. 349
  !!
  !! @par Revision History
  !!    A. J. Simmons, ECMWF, November 1981, original source
  !!    L. Kornblueh, MPI, May 1998, f90 rewrite
  !!    U. Schulzweida, MPI, May 1998, f90 rewrite
  !!
  SUBROUTINE auxhyb_chem_init(coord, ph,kdim,klen,nlev_in,                   &
                     pdelp,prdelp,plnph,plnpr,palpha )
    !<
    ! SUBROUTINE auxhyb                   
    ! Part of Module: mo_chem_init_coord
    ! based on mo_vertical_coord_table
    !
    ! Author: Jennifer Schroeter, KIT
    ! Initial Release: 2015-11-15                
    ! Modifications:
    !>
    TYPE(t_art_chem_init_coord), INTENT(in) :: &
     & coord
    INTEGER ,INTENT(in)  :: kdim, klen, nlev_in
    REAL(wp),INTENT(in)  :: ph(kdim, nlev_in+1)

    REAL(wp), INTENT(inout) :: pdelp (kdim, nlev_in  ), prdelp(kdim, nlev_in)
    REAL(wp), INTENT(inout) :: plnph (kdim, nlev_in+1), plnpr (kdim, nlev_in)
    REAL(wp), INTENT(inout) :: palpha(kdim, nlev_in  )

    REAL(wp) :: za, zd, zl, zr
    INTEGER  :: jk, jl, nlev, nlevp1

    nlev = nlev_in
    nlevp1 = nlev+1

    !-----
    ! Set pressure-level values or other top-level values

    DO jk = 1, coord%nplev
      zd = coord%delpr_chem_init(jk)
      zr = coord%rdelpr_chem_init(jk)
      zl = coord%rlnpr_chem_init(jk)
      za = coord%ralpha_chem_init(jk)
      DO jl = 1, klen
        pdelp(jl,jk) = zd
        prdelp(jl,jk) = zr
        plnpr(jl,jk) = zl
        palpha(jl,jk) = za
      END DO
    END DO

    ! Calculate hybrid-level values

    DO jk = coord%nplvp1, coord%nlmsgl

      IF (jk == 1 .AND. coord%vct_chem_init(1) == 0.0_wp ) THEN

        DO jl = 1, klen
          pdelp(jl,jk)  = ph(jl,jk+1) - ph(jl,jk)
          prdelp(jl,jk) = 1._wp/pdelp(jl,jk)
          palpha(jl,jk) = rd*LOG(2._wp)
          plnpr(jl,jk)  = 2._wp*palpha(jl,jk)
        END DO

      ELSE
        DO jl = 1, klen
          pdelp(jl,jk)  = ph(jl,jk+1) - ph(jl,jk)
          prdelp(jl,jk) = 1._wp/pdelp(jl,jk)
          plnpr(jl,jk)  = rd*LOG(ph(jl,jk+1)/ph(jl,jk))
          palpha(jl,jk) = rd - ph(jl,jk)*plnpr(jl,jk)*prdelp(jl,jk)
        END DO

      ENDIF
    END DO

    ! Set sigma-level values
    DO jk = coord%nlmsla, nlev
      zl = coord%rlnpr_chem_init(jk)
      za = coord%ralpha_chem_init(jk)
      DO jl = 1, klen
        pdelp(jl,jk)  = ph(jl,jk+1) - ph(jl,jk)
        prdelp(jl,jk) = 1._wp/pdelp(jl,jk)
        plnpr(jl,jk)  = zl
        palpha(jl,jk) = za
        
      END DO
    END DO

    DO jk = 2,nlevp1
      DO jl = 1, klen
        plnph(jl,jk) = LOG(ph(jl,jk))
      END DO
    END DO
  DO jl = 1, klen
        DO jk = 1, nlev_in
    ENDDO
    ENDDO

  END SUBROUTINE auxhyb_chem_init


  !-------------------------------------------------------------------------
  !>
  !! Calculates full- and half-level geopotential
  !!
  !! Method: Integrate the hydrostatic equation in the vertical
  !! to obtain full- or half-level values of geopotential.
  !!
  !! *geopot* is called during the calculation of adiabatic
  !! tendencies, prior to the calculation of the physical paramet-
  !! erizations, and in the post-processing.
  !!
  !! Parameters are
  !!  *pgeop_sfc*   *surface geopotential.
  !!  *pgeop_m*     *computed geopotential on full-levels
  !!  *pgeop_i*     *computed geopotential on half-levels
  !!  *ptv*         *virtual temperature.
  !!  *plnpr*       *logarithm of ratio of pressures, computed by *auxhyb*.
  !!  *palpha*      *for full-level values use *alpha* as computed by *auxhyb*.
  !!  *kdim*        *first dimension of 2-d arrays *phi,*
  !!                 *ptv,* *plnpr,* and *palpha.*
  !!  *kstart,kend* *start and end points for which calculation is performed.
  !!
  !! Required constants are obtained from module *mo_hyb*.
  !! The latter should have been initialized by a call of subroutine *inihyb*.
  !!
  !! Results are computed for *klen* consecutive points for *nlev* levels.
  !!
  !! The choice of full- or half-level values is determined
  !! by the specification of the input array *alpha.*
  !!
  !! External documentation of the model equations and the
  !! organization of the vertical calculation.
  !!
  !! @par Revision History
  !!    A. J. Simmons, ECMWF, January 1982, original source
  !!    L. Kornblueh, MPI, May 1998, f90 rewrite
  !!    U. Schulzweida, MPI, May 1998, f90 rewrite
  !!    H. Wan, MPI-M, 2006-08-10, included in m_vertical
  !!    H. Wan, MPI-M, 2006-08-15, half-level geopotential added to the output
  !!    H. Wan, MPI-M, 2007-07-19, calculation of the full-level geopotential
  !!                               simplified.
  !!    G. Zaengl, DWD, 2009-07-06, replace klen with kstart/kend for correct
  !!                                execution on nested domains
  !!    A. Seifert, DWD, 2010-06-21, add missing lower boundary by restructuring
  !!                                 the loop for half and full levels

  SUBROUTINE geopot_chem_init(coord, ptv,plnpr,palpha,pgeop_sfc,kdim,kstart,kend,nlev_in, &
                     pgeop_m, pgeop_i )
    !<
    ! SUBROUTINE geopot_chem_init                   
    ! Part of Module: mo_chem_init_coord
    ! based on mo_vertical_coord_table
    !
    ! Author: Jennifer Schroeter, KIT
    ! Initial Release: 2015-11-15                
    ! Modifications:
    !>
    TYPE(t_art_chem_init_coord) :: &
      &  coord
    INTEGER ,INTENT(in) :: kdim, kstart, kend, nlev_in

    REAL(wp) ,INTENT(in)    :: ptv   (kdim,nlev_in),     plnpr(kdim,nlev_in)
    REAL(wp) ,INTENT(in)    :: palpha(kdim,nlev_in), pgeop_sfc(kdim)
    REAL(wp) ,INTENT(inout) :: pgeop_m(kdim,nlev_in),   pgeop_i(kdim,nlev_in+1)

    INTEGER :: jk, jl, jkp, nlev, nlevp1

    nlev = nlev_in
    nlevp1 = nlev+1

    ! Integrate hydrostatic equation

    DO jl = kstart, kend
      pgeop_i(jl,nlevp1) = pgeop_sfc(jl)
      pgeop_m(jl,nlev)   = palpha(jl,nlev)*ptv(jl,nlev) + pgeop_sfc(jl)
    END DO

    ! half levels
    DO jk = nlev, 1, -1
      jkp = jk + 1
      DO jl = kstart, kend
        pgeop_i(jl,jk) = plnpr(jl,jk)*ptv(jl,jk) + pgeop_i(jl,jkp)
      END DO
    END DO

    ! full levels
    DO jk = coord%nlevm1, 1, -1
      jkp = jk + 1
      DO jl = kstart, kend
        pgeop_m(jl,jk)  = pgeop_i(jl,jkp) + palpha(jl,jk)*ptv(jl,jk)
      END DO
    END DO

  END SUBROUTINE geopot_chem_init
  
END MODULE mo_art_chem_init_coord


