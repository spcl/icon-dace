!
! Subroutine ml_ocean computes the surface temperature of a mixed
! layer slab ocean based only on the vertical energy budget.  currently
! there is no land/lakes, no ice, and no flux heat correction
!
! original purpose: use for a Radiative Convective Equilibrium APE exp
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

MODULE mo_ml_ocean

  USE mo_kind,               ONLY: wp
  USE mo_physical_constants, ONLY: rho_ref  

  IMPLICIT NONE
  PUBLIC :: ml_ocean

CONTAINS

  SUBROUTINE ml_ocean ( kproma, start_column, end_column, pdtime,pahflw,pahfsw,ptrflw,psoflw,ptsw )

    !
    ! Arguments
    !
    INTEGER, INTENT(IN)    :: kproma
    INTEGER, INTENT(IN)    :: start_column, end_column
    REAL(wp),INTENT(IN)    :: pdtime                   !< time step
    REAL(wp),INTENT(IN)    :: pahflw(kproma)           ! latent heat flux over water
    REAL(wp),INTENT(IN)    :: pahfsw(kproma)           ! sensible heat flux over water
    REAL(wp),INTENT(IN)    :: ptrflw(kproma)           ! LW flux over water
    REAL(wp),INTENT(IN)    :: psoflw(kproma)           ! SW flux over water
    REAL(wp),INTENT(INOUT) :: ptsw(kproma)             ! surface temperature of water
    !
    ! Local variables
    !
    INTEGER  :: jl
    REAL(wp) :: dmixsea,cpsea,zmixcap,zfluxw

    !$ACC DATA PRESENT(pahflw, pahfsw, ptrflw, psoflw, ptsw)

    !
    ! Parameters/constants set in ml_flux(ECHAM) originally:
    !  zdmix=50._dp
    !  zcpwater=3994._dp
    !  zrho_sea=1025._dp
    !  zmixcap=zrho_sea*zcpwater*zdmix
    !  zmonlen=2592000._dp
    !
    dmixsea   = 50._wp   ! [m] depth of mixed layer (zdmix)
    cpsea     = 3994._wp ! [J/(kg K)] specific heat capacity of water (what temp? looks low to me)

    zmixcap   = rho_ref*cpsea*dmixsea ! mixed layer heat capacity? 

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR PRIVATE(zfluxw)
    DO jl = start_column,  end_column

      zfluxw             = pahflw(jl) + pahfsw(jl) + ptrflw(jl) + psoflw(jl)
      ptsw(jl)           = ptsw(jl) + pdtime * zfluxw/zmixcap

    END DO
    !$ACC END PARALLEL

    !$ACC END DATA

  END SUBROUTINE ml_ocean

END MODULE mo_ml_ocean
