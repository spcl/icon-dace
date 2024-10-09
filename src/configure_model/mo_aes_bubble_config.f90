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

! @brief Module containing subroutines for the configuration of aes_bubble test case
!
! @literature cr2021_08_03)jsr for further documentation

MODULE mo_aes_bubble_config
                                                                                     
  USE mo_exception            ,ONLY: message, message_text, print_value, warning, finish                                                                                 
  USE mo_kind                 ,ONLY: wp                                              
  USE mo_parallel_config      ,ONLY: nproma                                          
                                                                                     
  IMPLICIT NONE                                                                      
                                                                                     
  PRIVATE                                                                            
                                                                                     
  ! configuration                                                                    
  PUBLIC ::         aes_bubble_config   !< user specified configuration parameters
  PUBLIC ::    init_aes_bubble_config   !< allocate and initialize aes_bubble_config
  PUBLIC ::    eval_aes_bubble_config   !< evaluate aes_bubble_config
  PUBLIC ::   print_aes_bubble_config   !< print out

  TYPE t_aes_bubble_config
     ! For a detailed description, see cr2021_08_03_jsr
     REAL(wp) :: psfc      !< surface temperature
     REAL(wp) :: t0        !< temperature at bottom of atmosphere (maybe different from surface temp)
     REAL(wp) :: t_am      !< absolute minimum of temperature (default: 180K)
     REAL(wp) :: gamma0    !< lapse rate in lowest atmospheric part
     REAL(wp) :: z0        !< below z0 lapse rate gamma0 is applied, above gamma1
     REAL(wp) :: gamma1    !< lapse rate above z0. However, temperature cannot fall below t_am
     REAL(wp) :: t_perturb !< maximum temperature perturbation in center of Gaussians
     REAL(wp) :: relhum_bg !< background relative humidity
     REAL(wp) :: relhum_mx !< maximum relative humidity
     REAL(wp) :: hw_x      ! half width in x-direction in meters
     REAL(wp) :: hw_z      ! half width in z-direction in meters
     REAL(wp) :: x_center  ! deviation of maximum of Gaussian in x direction
     LOGICAL  :: lgaussxy  ! .TRUE., if half width calculated for x-direction and x_center is applied also to y direction
  END TYPE t_aes_bubble_config

  TYPE(t_aes_bubble_config), TARGET :: aes_bubble_config

CONTAINS

  SUBROUTINE init_aes_bubble_config
    aes_bubble_config%psfc       = 101325._wp
    aes_bubble_config%t_am       = 180._wp
    aes_bubble_config%t0         = 303.15_wp
    aes_bubble_config%gamma0     = 0.009_wp
    aes_bubble_config%z0         = 3000._wp
    aes_bubble_config%gamma1     = 0.00001_wp
    aes_bubble_config%t_perturb  = 3._wp
    aes_bubble_config%relhum_bg  = 0.7_wp
    aes_bubble_config%relhum_mx  = 0.95_wp
    aes_bubble_config%hw_x       = 12500._wp
    aes_bubble_config%hw_z       = 500._wp
    aes_bubble_config%x_center   = 0._wp
    aes_bubble_config%lgaussxy   = .FALSE.
  END SUBROUTINE init_aes_bubble_config

  SUBROUTINE eval_aes_bubble_config
    CHARACTER(LEN=*), PARAMETER :: routine='eval_aes_bubble_config'
    REAL(wp), POINTER :: t_am,        psfc,      t0,           gamma0,      &
                       & gamma1,      z0,        t_perturb,    relhum_bg,   &
                       & relhum_mx,   hw_x,      hw_z,         x_center

    t_am        => aes_bubble_config%t_am
    psfc        => aes_bubble_config%psfc
    t0          => aes_bubble_config%t0
    gamma0      => aes_bubble_config%gamma0
    gamma1      => aes_bubble_config%gamma1
    z0          => aes_bubble_config%z0
    t_perturb   => aes_bubble_config%t_perturb
    relhum_bg   => aes_bubble_config%relhum_bg
    relhum_mx   => aes_bubble_config%relhum_mx
    hw_x        => aes_bubble_config%hw_x
    hw_z        => aes_bubble_config%hw_z
    x_center    => aes_bubble_config%x_center

    CALL message    ('','')                                                          
    CALL message    ('','------------------------------------------------------------------------')
    CALL message    ('','')                                                          
    CALL message    ('','Effective input to the aes bubble experiment on a torus')
    CALL message    ('','=======================================================')
    CALL message    ('','')                                                          
 
    IF (t_am <= 0._wp) THEN
       WRITE (message_text,*) 'Absolute minimum temperature in atmosphere t_am=',t_am, &
            & ' but cannot be below or equal to 0K'
       CALL finish(routine,message_text)
    ELSE
       WRITE (message_text,*) 'Absolute minimum temperature in atmosphere t_am=',t_am
       CALL message('',message_text)
    END IF
    ! t0 can be arbitrary, because the formula sets it to values >= t_am
    IF (psfc <= 0._wp) THEN
       WRITE (message_text,*) 'Surface pressure psfc=',psfc,' but cannot be below or equal to zero'
       CALL finish(routine,message_text)
    ELSE
       WRITE (message_text,*) 'Surface pressure psfc=',psfc
       CALL message('',message_text)
    END IF
    IF (gamma0 < 0._wp) THEN
       WRITE (message_text,*) 'Lapse rate gamma0=',gamma0,' but cannot be below 0'
       CALL finish(routine,message_text)
    ELSE
       WRITE (message_text,*) 'Lapse rate gamm0=',gamma0
       CALL message('',message_text)
    END IF
    IF (gamma1 < 0._wp) THEN
       WRITE (message_text,*) 'Lapse rate gamma1=',gamma1,' but cannot be below 0'
       CALL finish(routine,message_text)
    ELSE
       WRITE (message_text,*) 'Lapse rate gamm0=',gamma1
       CALL message('',message_text)
    END IF
    IF (z0 <= 0._wp) THEN
       WRITE (message_text,*) 'Altitude z0=',z0,' but cannot be below 0'
       CALL finish(routine,message_text)
    ELSE
       WRITE (message_text,*) 'Altitude z0=',z0
       CALL message('',message_text)
    END IF
    IF (t_perturb < 0._wp) THEN
       WRITE (message_text,*) 'Temperature perturbation in atmosphere t_perturb=',t_perturb, &
            & ' but cannot be below 0K'
       CALL finish(routine,message_text)
    ELSE
       WRITE (message_text,*) 'Temperature perturbation in atmosphere t_perturb=',t_perturb
       CALL message('',message_text)
    END IF
    IF (relhum_bg < 0._wp .or. relhum_bg > 1._wp) THEN
       WRITE (message_text,*) 'Background relative humidity relhum_bg=',relhum_bg, &
            ' but must be in [0,1]'
       CALL finish(routine,message_text)
    ELSE
       WRITE (message_text,*) 'Background relative humidity relhum_bg=',relhum_bg
       CALL message('',message_text)
    END IF
    IF (relhum_mx < 0._wp .or. relhum_mx > 1._wp) THEN
       WRITE (message_text,*) 'Relative humidity at maximum of Gaussians relhum_mx=', &
            & relhum_mx,' but must be in [0,1]'
       CALL finish(routine,message_text)
    ELSE
       WRITE (message_text,*) 'Relative humidity at maximum of Gaussians relhum_mx=',relhum_mx
       CALL message('',message_text)
    END IF
    IF (hw_x < 0._wp) THEN
       WRITE (message_text,*) 'Half width in horizontal direction hw_x=', hw_x,' but cannot be below 0'
       CALL finish(routine,message_text)
    ELSE
       WRITE (message_text,*) 'Half width in horizontal direction hw_x=', hw_x
       CALL message('',message_text)
    END IF
    IF (hw_z < 0._wp) THEN
       WRITE (message_text,*) 'Half width in vertical direction hw_z=',hw_z,' but cannot be below 0'
       CALL finish(routine,message_text)
    ELSE
       WRITE (message_text,*) 'Half width in vertical direction hw_z=',hw_z
       CALL message('',message_text)
    END IF
    CALL message ('','')
  END SUBROUTINE eval_aes_bubble_config

!!!==============================================================================================

  SUBROUTINE print_aes_bubble_config
    CALL message    ('','========================================================================')                                                                      
    CALL message    ('','')                                                          
    CALL message    ('','Aes bubble experiment configuration (on a Torus)')
    CALL message    ('','================================================')
    CALL message    ('','')                                                          

    CALL print_value ('    aes_bubble_config%t_am       ', aes_bubble_config%t_am      )
    CALL print_value ('    aes_bubble_config%psfc       ', aes_bubble_config%psfc      )
    CALL print_value ('    aes_bubble_config%t0         ', aes_bubble_config%t0        )
    CALL print_value ('    aes_bubble_config%gamma0     ', aes_bubble_config%gamma0    )
    CALL print_value ('    aes_bubble_config%gamma1     ', aes_bubble_config%gamma1    )
    CALL print_value ('    aes_bubble_config%z0         ', aes_bubble_config%z0        )
    CALL print_value ('    aes_bubble_config%t_perturb  ', aes_bubble_config%t_perturb )
    CALL print_value ('    aes_bubble_config%relhum_bg  ', aes_bubble_config%relhum_bg )
    CALL print_value ('    aes_bubble_config%relhum_mx  ', aes_bubble_config%relhum_mx )
    CALL print_value ('    aes_bubble_config%hw_x       ', aes_bubble_config%hw_x      )
    CALL print_value ('    aes_bubble_config%hw_z       ', aes_bubble_config%hw_z      )
    CALL print_value ('    aes_bubble_config%x_center   ', aes_bubble_config%x_center  )
    CALL print_value ('    aes_bubble_config%lgaussxz   ', aes_bubble_config%lgaussxy  )
    
  END SUBROUTINE print_aes_bubble_config

!!!==============================================================================================

END MODULE mo_aes_bubble_config
