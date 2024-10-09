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

! Data types and variables used by the rte-rrtmgp radiation package.
!
! This module contains
!  definition of data types for organising rte-rrtmgp diagnostics,
!  the actual variables that are declared of these types, and
!  subroutines for (de-)allocating memory for the variables.

MODULE mo_radiation_forcing
  USE mo_kind,                ONLY: wp
  USE mo_radiation_forcing_memory, ONLY: psf=>prm_radiation_forcing
  USE mo_physical_constants,  ONLY: stbo
  USE mo_radiation_config,    ONLY: lradforcing
  USE mtime,                  ONLY: no_of_sec_in_a_day

  IMPLICIT NONE
  PRIVATE

  PUBLIC          :: prepare_radiation_forcing, calculate_radiation_forcing

  CONTAINS

  !-----------------------------------------------------------------------------
  ! prepare_forcing: calculate intermediate quantities for forcing calculation
  !-----------------------------------------------------------------------------

  SUBROUTINE prepare_radiation_forcing(jg                                ,&
       & kproma           ,kbdim             ,klevp1           ,krow           ,&
       & pflx_dnlw        ,pflx_dnsw         ,pflx_dnlw_clr    ,pflx_dnsw_clr   )

    INTEGER, INTENT(in)     :: jg,     & !< domain index 
                             & kproma, & !< block length
                             & kbdim,  & !< maximum block length
                             & klevp1, & !< number of layers + 1
                             & krow      !< block index
    REAL(wp), INTENT(in)    :: &
         pflx_dnlw(kbdim,klevp1),    & !< Net dwnwrd LW flux [Wm2]
         pflx_dnsw(kbdim,klevp1),    & !< Net dwnwrd SW flux [Wm2] for forcing
         pflx_dnlw_clr(kbdim,klevp1),& !< Net dn LW flux (clear sky) [Wm2]
         pflx_dnsw_clr(kbdim,klevp1)   !< Net dn SW flux (clear sky) [Wm2]

    psf(jg)%emter_for(1:kproma,1:klevp1,krow)=pflx_dnlw(1:kproma,1:klevp1)
    psf(jg)%emtef_for(1:kproma,1:klevp1,krow)=pflx_dnlw_clr(1:kproma,1:klevp1)
    psf(jg)%trsol_for(1:kproma,1:klevp1,krow)=pflx_dnsw(1:kproma,1:klevp1)
    psf(jg)%trsof_for(1:kproma,1:klevp1,krow)=pflx_dnsw_clr(1:kproma,1:klevp1)

  END SUBROUTINE prepare_radiation_forcing
  !-----------------------------------------------------------------------------
  ! calculate_forcing: calculate radiative forcing 
  !-----------------------------------------------------------------------------

  SUBROUTINE calculate_radiation_forcing( jg                 &
       & ,jcs                ,jce              ,kbdim              &
       & ,klevp1             ,krow             ,pi0                &
       & ,pconvfact          ,pflxs            ,pflxs0             &
       & ,pflxt              ,pflxt0           ,pemiss             &
       & ,ptsfctrad          ,pztsnew                              )
    INTEGER, INTENT(in)          :: jg,     & !< domain index 
                                  & jcs,jce, & !< block length
                                  & kbdim,  & !< maximum block length
                                  & klevp1, & !< number of layers + 1
                                  & krow      !< block index
    REAL(wp), INTENT(in)         :: &
         & pi0(kbdim),                 & !< solar irradiance at top of atmosphere
         & pconvfact(kbdim,klevp1-1),  & !< conversion factor radiation flux -> heating rate
         & pflxs(kbdim,klevp1),        & !< short wave net radiation flux (all sky)
         & pflxs0(kbdim,klevp1),       & !< short wave net radiation flux (clear sky)
         & pflxt(kbdim,klevp1),        & !< long wave net radiation flux (all sky)
         & pflxt0(kbdim,klevp1),       & !< long wave net radiation flux (clear sky)
         & pemiss(kbdim),              & !< emissivity at surface
         & ptsfctrad(kbdim),           & !< surface temperature at radiation time step
         & pztsnew(kbdim)                !< surface temperature

    INTEGER                     :: jk, jl, klev
    REAL(wp)                    :: zemter_for(kbdim,klevp1), zemtef_for(kbdim,klevp1), &
                                 & ztrsol_for(kbdim,klevp1), ztrsof_for(kbdim,klevp1)
    REAL(wp)                    :: zflxs_all_for(kbdim,klevp1), &
                                 & zflxs_clear_for(kbdim,klevp1), &
                                 & zflxt_all_for(kbdim,klevp1), &
                                 & zflxt_clear_for(kbdim,klevp1)
    REAL(wp)                    :: dlwem_o_dtg(kbdim)
    REAL(wp)                    :: zdtdt_sw, zdtdt_all_for_sw, &
         zdtdt_lw, zdtdt_all_for_lw

    klev=klevp1-1
    zemter_for(jcs:jce,1:klevp1)=psf(jg)%emter_for(jcs:jce,1:klevp1,krow)
    zemtef_for(jcs:jce,1:klevp1)=psf(jg)%emtef_for(jcs:jce,1:klevp1,krow)
    ztrsol_for(jcs:jce,1:klevp1)=psf(jg)%trsol_for(jcs:jce,1:klevp1,krow)
    ztrsof_for(jcs:jce,1:klevp1)=psf(jg)%trsof_for(jcs:jce,1:klevp1,krow)

    IF (lradforcing(1)) THEN
      !--- Radiation flux forcing:
      zflxs_all_for(jcs:jce,1:klevp1)=SPREAD(pi0(jcs:jce),2,klevp1)* &
           ztrsol_for(jcs:jce,1:klevp1)
      zflxs_clear_for(jcs:jce,1:klevp1)=SPREAD(pi0(jcs:jce),2,klevp1)* &
           ztrsof_for(jcs:jce,1:klevp1)
      !  forcing solar wave length bands:
      DO jk = 1, klev
        DO jl = jcs, jce
          zdtdt_sw=pconvfact(jl,jk)*(pflxs(jl,jk+1)-pflxs(jl,jk))
          zdtdt_all_for_sw=pconvfact(jl,jk)*(zflxs_all_for(jl,jk+1)-zflxs_all_for(jl,jk))
          psf(jg)%netht_sw(jl,jk,krow) = REAL(no_of_sec_in_a_day,wp)*(zdtdt_sw-zdtdt_all_for_sw)
        ENDDO
      END DO
      DO jk = 1, klevp1
        psf(jg)%d_aflx_sw(jcs:jce,jk,krow)  = pflxs(jcs:jce,jk) - zflxs_all_for(jcs:jce,jk)
        psf(jg)%d_aflx_swc(jcs:jce,jk,krow) = pflxs0(jcs:jce,jk)- zflxs_clear_for(jcs:jce,jk)
      END DO
      psf(jg)%fsw_total_top(jcs:jce,krow)=psf(jg)%d_aflx_sw(jcs:jce,1,krow)
      psf(jg)%fsw_total_sur(jcs:jce,krow)=psf(jg)%d_aflx_sw(jcs:jce,klevp1,krow)
      psf(jg)%fsw_clear_top(jcs:jce,krow)=psf(jg)%d_aflx_swc(jcs:jce,1,krow)
      psf(jg)%fsw_clear_sur(jcs:jce,krow)=psf(jg)%d_aflx_swc(jcs:jce,klevp1,krow)
    END IF

    IF (lradforcing(2)) THEN
      !--- Radiation flux forcing:
      zflxt_all_for(jcs:jce,1:klev)=zemter_for(jcs:jce,1:klev)
      ! 
      ! in the following formulae, the order of calculation has to be as is
      ! if not, for an aerosol free atmosphere the forcing is not exactly 0
      ! because of numeric effects (see the corresponding formulae in radheat)
      ! NEW adjustment of radiation at surface by Tayler expansion
      dlwem_o_dtg(jcs:jce) =  pemiss(jcs:jce)*4._wp*stbo*pztsnew(jcs:jce)**3
!!$      zflxt_all_for(jcs:jce,klevp1)   = zemter_for(jcs:jce,klevp1)+ &
!!$           pemiss(jcs:jce)*stbo*ptsfctrad(jcs:jce)**4
!!$      zflxt_all_for(jcs:jce,klevp1)   = zflxt_all_for(jcs:jce,klevp1)- &
!!$           pemiss(jcs:jce)*stbo*pztsnew(jcs:jce)**4
      zflxt_all_for(jcs:jce,klevp1)   = zemter_for(jcs:jce,klevp1) &
           & - dlwem_o_dtg(jcs:jce) * (pztsnew(jcs:jce) - ptsfctrad(jcs:jce))
      zflxt_clear_for(jcs:jce,1:klev) = zemtef_for(jcs:jce,1:klev)
!!$      zflxt_clear_for(jcs:jce,klevp1) = zemtef_for(jcs:jce,klevp1) + &
!!$           pemiss(jcs:jce)*stbo*ptsfctrad(jcs:jce)**4
!!$      zflxt_clear_for(jcs:jce,klevp1) = zflxt_clear_for(jcs:jce,klevp1)- &
!!$           pemiss(jcs:jce)*stbo*pztsnew(jcs:jce)**4
      zflxt_clear_for(jcs:jce,klevp1) = zemtef_for(jcs:jce,klevp1) &
           & - dlwem_o_dtg(jcs:jce) * (pztsnew(jcs:jce) - ptsfctrad(jcs:jce))
      DO jk = 1,klev
        DO jl = jcs, jce
          zdtdt_lw=pconvfact(jl,jk)*(pflxt(jl,jk+1)-pflxt(jl,jk))
          zdtdt_all_for_lw=pconvfact(jl,jk)*(zflxt_all_for(jl,jk+1)-zflxt_all_for(jl,jk))
          psf(jg)%netht_lw(jl,jk,krow) = REAL(no_of_sec_in_a_day,wp)*(zdtdt_lw-zdtdt_all_for_lw)
        END DO
      END DO
      DO jk = 1, klevp1
        psf(jg)%d_aflx_lw(jcs:jce,jk,krow)  = pflxt(jcs:jce,jk)-zflxt_all_for(jcs:jce,jk)
        psf(jg)%d_aflx_lwc(jcs:jce,jk,krow) = pflxt0(jcs:jce,jk)-zflxt_clear_for(jcs:jce,jk)
      END DO
      psf(jg)%flw_total_top(jcs:jce,krow) = psf(jg)%d_aflx_lw(jcs:jce,1,krow)
      psf(jg)%flw_total_sur(jcs:jce,krow) = psf(jg)%d_aflx_lw(jcs:jce,klevp1,krow)
      psf(jg)%flw_clear_top(jcs:jce,krow) = psf(jg)%d_aflx_lwc(jcs:jce,1,krow)
      psf(jg)%flw_clear_sur(jcs:jce,krow) = psf(jg)%d_aflx_lwc(jcs:jce,klevp1,krow)
    END IF
  END SUBROUTINE calculate_radiation_forcing

END MODULE mo_radiation_forcing
