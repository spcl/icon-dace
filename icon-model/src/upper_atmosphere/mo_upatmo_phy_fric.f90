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

! Computation of frictional heating in a physically meaningful way.

MODULE mo_upatmo_phy_fric

  USE mo_kind,               ONLY: wp
  USE mo_impl_constants,     ONLY: SUCCESS
  USE mo_physical_constants, ONLY: argas, tmelt

  IMPLICIT NONE
  
  PUBLIC :: fric_heat

CONTAINS
  
  !>
  !! Compute frictional heating
  !!
  !! Literature:
  !! - Gill, A. E. (1982) Atmosphere-ocean dynmaics. Academic Press, London, 4 edn.
  !!
  SUBROUTINE fric_heat(jcs, jce, kbdim, klev, ptm1, ptvm1, pum1, pvm1, papm1, paphm1, grav, amu, cp, ptte_fc, &
    &                  opt_istartlev, opt_iendlev, opt_ldiss_from_heatdiff)

    ! in/out variables
    INTEGER,  INTENT(IN)  :: jcs, jce, kbdim, klev
    REAL(wp), INTENT(IN)  :: ptm1(kbdim,klev)       ! temperature           [T]
    REAL(wp), INTENT(IN)  :: ptvm1(kbdim,klev)      ! virtual temperature   [TV]
    REAL(wp), INTENT(IN)  :: pum1(kbdim,klev)       ! u-velocity            [u]
    REAL(wp), INTENT(IN)  :: pvm1(kbdim,klev)       ! v-velocity            [v]
    REAL(wp), INTENT(IN)  :: papm1(kbdim,klev)      ! full-level pressure   [pf]
    REAL(wp), INTENT(IN)  :: paphm1(kbdim,klev+1)   ! half-level pressure   [ph]
    REAL(wp), INTENT(IN)  :: grav(kbdim,klev)       ! gravity acceleration  [g]
    REAL(wp), INTENT(IN)  :: amu(kbdim,klev)        ! molar weight of air   [am]
    REAL(wp), INTENT(IN)  :: cp(kbdim,klev)         ! specific heat of air  [cp]

    REAL(wp), INTENT(OUT) :: ptte_fc(kbdim,klev)    ! temperature tendency due to frictional heating [dq/dt]

    INTEGER,  OPTIONAL, INTENT(IN) :: opt_istartlev, opt_iendlev ! optional vertical start and end indices
    LOGICAL,  OPTIONAL, INTENT(IN) :: opt_ldiss_from_heatdiff    ! heat source from heat diffusion

    ! local variables
    INTEGER  :: jl, jk, istartlev, iendlev, iendlevm1, istartlevp1
    REAL(wp) :: zrstar, zrdz, zdudz, zdvdz, zcoef, ztempdz
    REAL(wp) :: zth, zgvh, zmah, ztvh, zrhoh
    REAL(wp), ALLOCATABLE :: zgmurhoh(:,:)
    LOGICAL  :: ldiss_from_heatdiff

    REAL(wp), PARAMETER :: inv_tmelt = 1._wp / tmelt   ! (tmelt=273.15_wp)
    REAL(wp), PARAMETER :: inv_Pr    = 1._wp / 0.72_wp ! inverse Prandtl number

    !---------------------------------------------------------

    ! please do not limit range of assignment 
    ! (e.g., ptte_fc(jcs:jce,istartlev:iendlev) = 0._wp)), 
    ! since tendencies have attribute INTENT(OUT)
    ptte_fc(:,:) = 0._wp

    ! determine start and end indices of vertical grid layers,
    ! for which tendencies should be computed
    IF (PRESENT(opt_istartlev)) THEN
      istartlev = MIN(MAX(1, opt_istartlev), klev)
    ELSE
      istartlev = 1
    ENDIF

    IF (PRESENT(opt_iendlev)) THEN
      iendlev = MIN(MAX(1, opt_iendlev), klev)
    ELSE
      iendlev = klev
    ENDIF

    IF (istartlev >= iendlev) RETURN 

    IF (PRESENT(opt_ldiss_from_heatdiff)) THEN
      ldiss_from_heatdiff = opt_ldiss_from_heatdiff
    ELSE
      ldiss_from_heatdiff = .FALSE.
    ENDIF
    
    istartlevp1 = istartlev + 1
    iendlevm1   = iendlev - 1

    zrstar = 1.e3_wp * argas  ! universal gas constant 8314 J/kmol/K

    ALLOCATE(zgmurhoh(kbdim,istartlevp1:iendlev))

    DO jk = istartlevp1, iendlev
      DO jl = jcs, jce
        zth   = 0.5_wp * ( ptm1(jl,jk) + ptm1(jl,jk-1) )
        zgvh  = 0.5_wp * ( grav(jl,jk) + grav(jl,jk-1) )
        zmah  = 0.5_wp * ( amu(jl,jk) + amu(jl,jk-1) )
        ztvh  = 0.5_wp * ( ptvm1(jl,jk) + ptvm1(jl,jk-1) )
        zrhoh = paphm1(jl,jk) * zmah / ( zrstar * ztvh )
        zgmurhoh(jl,jk) = 1.87E-5_wp * EXP( 0.69_wp * LOG( inv_tmelt * zth ) ) &
          &             * zgvh * zrhoh
      ENDDO  !jl
    ENDDO  !jk

    DO jk = istartlevp1, iendlevm1
      DO jl = jcs, jce
        zrdz  = 1._wp / ( papm1(jl,jk-1) - papm1(jl,jk+1) )
        zdudz = ( pum1(jl,jk-1) - pum1(jl,jk+1) ) * zrdz
        zdvdz = ( pvm1(jl,jk-1) - pvm1(jl,jk+1) ) * zrdz
        zcoef = 0.5_wp * ( zgmurhoh(jl,jk) + zgmurhoh(jl,jk+1) ) * grav(jl,jk) / cp(jl,jk)
        ptte_fc(jl,jk) = zcoef * ( zdudz * zdudz + zdvdz * zdvdz )
      ENDDO  !jl
    ENDDO  !jk

    DO jl = jcs, jce
      zrdz  = 1._wp / ( papm1(jl,istartlev) - papm1(jl,istartlevp1) )
      zdudz = ( pum1(jl,istartlev) - pum1(jl,istartlevp1) ) * zrdz 
      zdvdz = ( pvm1(jl,istartlev) - pvm1(jl,istartlevp1) ) * zrdz
      zcoef = zgmurhoh(jl,istartlevp1) * grav(jl,istartlev) / cp(jl,istartlev)
      ptte_fc(jl,istartlev) = zcoef * ( zdudz * zdudz + zdvdz * zdvdz )

      zrdz  = 1._wp / ( papm1(jl,iendlev) - papm1(jl,iendlevm1) )
      zdudz = ( pum1(jl,iendlev) - pum1(jl,iendlevm1) ) * zrdz
      zdvdz = ( pvm1(jl,iendlev) - pvm1(jl,iendlevm1) ) * zrdz
      zcoef = zgmurhoh(jl,iendlev) * grav(jl,iendlev) / cp(jl,iendlev)
      ptte_fc(jl,iendlev) = zcoef * ( zdudz * zdudz + zdvdz * zdvdz )
    ENDDO  !jl

    IF (ldiss_from_heatdiff) THEN

      ! compute heat source from heat diffusion

      DO jk = istartlevp1, iendlevm1
        DO jl = jcs, jce
          ztempdz = ( ptm1(jl,jk-1) - ptm1(jl,jk+1) ) / &
            &       ( papm1(jl,jk-1) - papm1(jl,jk+1) )
          zcoef   = 0.5_wp * inv_Pr * ( zgmurhoh(jl,jk) + zgmurhoh(jl,jk+1) ) * grav(jl,jk) / cp(jl,jk)
          ptte_fc(jl,jk) = ptte_fc(jl,jk) + zcoef * ztempdz**2
        ENDDO  !jl
      ENDDO  !jk
      
      DO jl = jcs, jce
        ztempdz = ( ptm1(jl,istartlev) - ptm1(jl,istartlevp1) ) / &
          &       ( papm1(jl,istartlev) - papm1(jl,istartlevp1) )
        zcoef   = inv_Pr * zgmurhoh(jl,istartlevp1) * grav(jl,istartlev) / cp(jl,istartlev)
        ptte_fc(jl,istartlev) = ptte_fc(jl,istartlev) + zcoef * ztempdz**2
        
        ztempdz = ( pum1(jl,iendlev) - pum1(jl,iendlevm1) ) / &
          &       ( papm1(jl,iendlev) - papm1(jl,iendlevm1) )
        zcoef   = inv_Pr * zgmurhoh(jl,iendlev) * grav(jl,iendlev) / cp(jl,iendlev)
        ptte_fc(jl,iendlev) = ptte_fc(jl,iendlev) + zcoef * ztempdz**2
      ENDDO  !jl

    ENDIF

    IF (ALLOCATED(zgmurhoh)) DEALLOCATE(zgmurhoh)

  END SUBROUTINE fric_heat

END MODULE mo_upatmo_phy_fric
