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

MODULE mo_diffusion_config

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, message_text, finish
  USE mo_impl_constants,      ONLY: max_dom

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: t_diffusion_config, diffusion_config  !< derived type and variable
  PUBLIC :: configure_diffusion                   !< subroutine

  !--------------------------------------------------------------------------
  ! Basic configuration setup for diffusion
  !--------------------------------------------------------------------------
  TYPE t_diffusion_config

    ! variables from namelist

    INTEGER :: hdiff_order  ! order of horizontal diffusion
                            ! -1: no diffusion
                            ! 2: 2nd order linear diffusion on all vertical levels
                            ! 4: 4th order linear diffusion on all vertical levels
                            ! 5: Smagorinsky diffusion with optional fourth-order background diffusion

    REAL(wp) :: hdiff_efdt_ratio      ! ratio of e-folding time to (2*)time step
    REAL(wp) :: hdiff_w_efdt_ratio    ! ratio of e-folding time to time step for w diffusion (NH only)
    REAL(wp) :: hdiff_min_efdt_ratio  ! minimum value of hdiff_efdt_ratio 
                                      ! (for upper sponge layer)
    REAL(wp) :: hdiff_tv_ratio        ! the ratio of diffusion coefficient: temp:mom
    REAL(wp) :: hdiff_smag_fac        ! scaling factor for Smagorinsky diffusion at height hdiff_smag_z and below
    REAL(wp) :: hdiff_smag_fac2       ! scaling factor for Smagorinsky diffusion at height hdiff_smag_z2
    REAL(wp) :: hdiff_smag_fac3       ! scaling factor for Smagorinsky diffusion at height hdiff_smag_z3
    REAL(wp) :: hdiff_smag_fac4       ! scaling factor for Smagorinsky diffusion at height hdiff_smag_z4 and above
    REAL(wp) :: hdiff_smag_z          ! height up to which hdiff_smag_fac is used, start of linear profile
    REAL(wp) :: hdiff_smag_z2         ! height of hdiff_smag_fac2, end of linear and start of quadratic profile
    REAL(wp) :: hdiff_smag_z3         ! height of hdiff_smag_fac3, to define quadratic profile
    REAL(wp) :: hdiff_smag_z4         ! height from which hdiff_smag_fac4, end of quadratic profile
    REAL(wp) :: hdiff_multfac         ! multiplication factor of normalized diffusion
                                      ! coefficient for nested domains
    INTEGER  :: itype_vn_diffu        ! options for discretizing the Smagorinsky momentum diffusion
    INTEGER  :: itype_t_diffu         ! options for discretizing the Smagorinsky temperature diffusion

    LOGICAL :: lhdiff_temp   ! if .TRUE., apply horizontal diffusion to temp.
    LOGICAL :: lhdiff_vn     ! if .TRUE., apply horizontal diffusion to momentum.
    LOGICAL :: lhdiff_w      ! if .TRUE., apply horizontal diffusion to vertical momentum.
    LOGICAL :: lsmag_3d      ! if .TRUE., compute 3D Smagorinsky diffusion coefficient.
    LOGICAL :: lhdiff_smag_w ! if .TRUE., apply additional Smagorinsky diffusion to vertical momentum.

    ! variables not from namelist

    REAL(wp) :: k6, k4, k2, k4w  ! numerical diffusion coefficients
                                 ! Values for these parameters are not directly
                                 ! specified by the user, but derived from the ratio 
                                 ! between the e-folding time and the model time step
                                 ! (hdiff_efdt_ratio above), and the horizontal 
                                 ! resolution of the model

  END TYPE t_diffusion_config
  !>
  !!
  TYPE(t_diffusion_config) :: diffusion_config(max_dom)

CONTAINS
  !>
  !!
  SUBROUTINE configure_diffusion( n_dom, parent_id )

    INTEGER, INTENT(IN) :: n_dom
    INTEGER, INTENT(IN) :: parent_id(:)

    INTEGER  :: jg, jgp
    REAL(wp) :: tmp

    CHARACTER(len=*), PARAMETER :: &
      routine = 'mo_diffusion_config:configure_diffusion'

    !-----------------------------------------------------------
    ! Compute diffusion coefficients
    !-----------------------------------------------------------

    IF (ANY(diffusion_config(1:n_dom)%hdiff_efdt_ratio<=0._wp)) THEN

      diffusion_config(:)%k2 = 0._wp
      diffusion_config(:)%k4 = 0._wp
      diffusion_config(:)%k6 = 0._wp

      CALL message(TRIM(routine),'Background linear diffusion is '//&
                                 'switched off in all domains')

    ELSE

      tmp = 1._wp/(diffusion_config(1)%hdiff_efdt_ratio*8._wp)
      diffusion_config(1)%k2 = tmp
      tmp = 1._wp/(diffusion_config(1)%hdiff_efdt_ratio*64._wp)
      diffusion_config(1)%k4 = tmp
      tmp = 1._wp/(diffusion_config(1)%hdiff_efdt_ratio*512._wp)
      diffusion_config(1)%k6 = tmp

      tmp = 1._wp/(diffusion_config(1)%hdiff_w_efdt_ratio*36._wp)
      diffusion_config(:)%k4w = tmp

      DO jg = 2, n_dom

         jgp = parent_id(jg)

         diffusion_config(jg)%k2 = diffusion_config(jgp)%k2 * diffusion_config(jg)%hdiff_multfac
         diffusion_config(jg)%k4 = diffusion_config(jgp)%k4 * diffusion_config(jg)%hdiff_multfac
         diffusion_config(jg)%k6 = diffusion_config(jgp)%k6 * diffusion_config(jg)%hdiff_multfac
      ENDDO

    ENDIF


    ! produce some log-output
    !
    DO jg =1,n_dom
      SELECT CASE( diffusion_config(jg)%hdiff_order )
      CASE(-1)
        WRITE(message_text,'(a,i2.2)') 'Horizontal diffusion '//&
                                       'switched off for domain ', jg
        CALL message(routine, message_text)

      CASE(2,4,5)
        ! do nothing

      CASE DEFAULT
        WRITE(message_text,'(a,i2.2,a)') 'Error: Invalid choice for hdiff_order '//&
          &                              'for domain ', jg,                        &
          &                              '. Choose from -1, 2, 4, and 5.'
        CALL finish(routine, message_text)
      END SELECT

      IF ( diffusion_config(jg)%hdiff_efdt_ratio<=0._wp ) THEN
        WRITE(message_text,'(a,i2.2)') 'No horizontal background diffusion is used '//&
                                       'for domain ', jg
        CALL message(routine, message_text)
      ENDIF
    ENDDO

  END SUBROUTINE configure_diffusion

END MODULE mo_diffusion_config
