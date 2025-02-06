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

! Configuration of the parameterization for vertical diffusion,
! that is used in the AES physics package.
!
! References:
!     Angevine, W. M., Jiang, H., & Mauritsen T. (2010).
!           Performance of an eddy diffusivity- mass flux scheme for shallow cumulus boundary layers.
!           Monthly Weather Review, 138(7), 2895-2912. https://doi.org/10.1175/2010MWR3142.1
!     Mauritsen, T., & Svensson, G. (2007).
!           Observations of stably stratified shear-driven atmospheric turbulence at low and high Richardson numbers.
!           Journal of the Atmospheric Sciences, 64(2), 645-655. https://doi.org/10.1175/JAS3856.1

MODULE mo_aes_vdf_config

  USE mo_turb_vdiff_config    ,ONLY: t_vdiff_config, vdiff_config_init, &
    vdiff_config_update, vdiff_config_check
  USE mo_impl_constants       ,ONLY: max_dom
  USE mo_grid_config          ,ONLY: n_dom
  USE mo_exception,            ONLY: message, print_value

  IMPLICIT NONE

  PRIVATE

  PUBLIC ::                     name   !< name for this unit

  ! configuration
  PUBLIC ::         aes_vdf_config     !< user specified configuration parameters
  PUBLIC ::    init_aes_vdf_config     !< allocate and initialize aes_vdf_config
  PUBLIC ::    eval_aes_vdf_config     !< evaluate aes_vdf_config
  PUBLIC ::   print_aes_vdf_config     !< print out

  !>
  !! Name of this unit
  !!
  CHARACTER(LEN=*), PARAMETER :: name = 'aes_vdf'

  !>
  !! Configuration state vectors, for multiple domains/grids.
  !!
  TYPE(t_vdiff_config), TARGET :: aes_vdf_config(max_dom)

CONTAINS

  !----

  !>
  !! Initialize the global configuration state vector
  !!
  SUBROUTINE init_aes_vdf_config
    CALL vdiff_config_init(aes_vdf_config(:))
  END SUBROUTINE init_aes_vdf_config

  !>
  !! Evaluate additional derived parameters
  !!
  SUBROUTINE eval_aes_vdf_config
    CALL vdiff_config_update(aes_vdf_config(:))
    CALL vdiff_config_check(aes_vdf_config(:))
  END SUBROUTINE eval_aes_vdf_config

  !----

  !>
  !! Print out the user controlled configuration state
  !!
  SUBROUTINE print_aes_vdf_config(ng)
    !
    INTEGER, INTENT(in) :: ng
    !
    INTEGER             :: jg
    CHARACTER(LEN=2)    :: cg
    !
    CALL message    ('','')
    CALL message    ('','========================================================================')
    CALL message    ('','')
    CALL message    ('','AES vertical diffusion configuration')
    CALL message    ('','======================================')
    CALL message    ('','')
    !
    DO jg = 1,ng
       !
       WRITE(cg,'(i0)') jg
       !
       CALL message    ('','For domain '//cg)
       CALL message    ('','------------')
       CALL message    ('','')
       CALL print_value('    aes_vdf_config('//TRIM(cg)//')% lsfc_mom_flux  ',aes_vdf_config(jg)% lsfc_mom_flux  )
       CALL print_value('    aes_vdf_config('//TRIM(cg)//')% lsfc_heat_flux ',aes_vdf_config(jg)% lsfc_heat_flux )
       CALL message    ('','')
       CALL print_value('    aes_vdf_config('//TRIM(cg)//')% pr0            ',aes_vdf_config(jg)% pr0            )
       CALL print_value('    aes_vdf_config('//TRIM(cg)//')% f_tau0         ',aes_vdf_config(jg)% f_tau0         )
       CALL print_value('    aes_vdf_config('//TRIM(cg)//')% f_theta0       ',aes_vdf_config(jg)% f_theta0       )
       CALL print_value('    aes_vdf_config('//TRIM(cg)//')% c_f            ',aes_vdf_config(jg)% c_f            )
       CALL print_value('    aes_vdf_config('//TRIM(cg)//')% c_n            ',aes_vdf_config(jg)% c_n            )
       CALL print_value('    aes_vdf_config('//TRIM(cg)//')% c_e            ',aes_vdf_config(jg)% c_e            )
       CALL print_value('    aes_vdf_config('//TRIM(cg)//')% wmc            ',aes_vdf_config(jg)% wmc            )
       CALL print_value('    aes_vdf_config('//TRIM(cg)//')% fsl            ',aes_vdf_config(jg)% fsl            )
       CALL print_value('    aes_vdf_config('//TRIM(cg)//')% fbl            ',aes_vdf_config(jg)% fbl            )
       CALL print_value('    aes_vdf_config('//TRIM(cg)//')% lmix_max       ',aes_vdf_config(jg)% lmix_max       )
       CALL print_value('    aes_vdf_config('//TRIM(cg)//')% z0m_min        ',aes_vdf_config(jg)% z0m_min        )
       CALL print_value('    aes_vdf_config('//TRIM(cg)//')% z0m_ice        ',aes_vdf_config(jg)% z0m_ice        )
       CALL print_value('    aes_vdf_config('//TRIM(cg)//')% z0m_oce        ',aes_vdf_config(jg)% z0m_oce        )
       CALL print_value('    aes_vdf_config('//TRIM(cg)//')% turb           ',aes_vdf_config(jg)% turb           )
       CALL print_value('    aes_vdf_config('//TRIM(cg)//')% use_tmx        ',aes_vdf_config(jg)% use_tmx        )
       IF (aes_vdf_config(jg)% use_tmx) THEN
        CALL print_value('    aes_vdf_config('//TRIM(cg)//')% solver_type    ',aes_vdf_config(jg)% solver_type   )
        CALL print_value('    aes_vdf_config('//TRIM(cg)//')% energy_type    ',aes_vdf_config(jg)% energy_type   )
        CALL print_value('    aes_vdf_config('//TRIM(cg)//')% dissipation_factor',aes_vdf_config(jg)% dissipation_factor)
        CALL print_value('    aes_vdf_config('//TRIM(cg)//')% use_louis     ',aes_vdf_config(jg)% use_louis      )
        CALL print_value('    aes_vdf_config('//TRIM(cg)//')% louis_constant_b',aes_vdf_config(jg)% louis_constant_b)
       END IF
       CALL print_value('    aes_vdf_config('//TRIM(cg)//')% smag_constant  ',aes_vdf_config(jg)% smag_constant  )
       CALL print_value('    aes_vdf_config('//TRIM(cg)//')% turb_prandtl   ',aes_vdf_config(jg)% turb_prandtl   )
       CALL print_value('    aes_vdf_config('//TRIM(cg)//')% rturb_prandtl  ',aes_vdf_config(jg)% rturb_prandtl  )
       CALL print_value('    aes_vdf_config('//TRIM(cg)//')% km_min         ',aes_vdf_config(jg)% km_min         )
       CALL print_value('    aes_vdf_config('//TRIM(cg)//')% max_turb_scale ',aes_vdf_config(jg)% max_turb_scale )
       CALL print_value('    aes_vdf_config('//TRIM(cg)//')% min_sfc_wind   ',aes_vdf_config(jg)% min_sfc_wind   )
       CALL message    ('','')
       !
    END DO
    !
  END SUBROUTINE print_aes_vdf_config

  !----

END MODULE mo_aes_vdf_config
