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

! mm: modularized aggregation scheme following Maerz et al. (2020)

MODULE mo_aggregates


     USE mo_kind,    ONLY : wp

     USE mo_bgc_constants,  ONLY : pi, g, rhoref_water

     USE mo_control_bgc,    ONLY : dtbgc, bgc_nproma, bgc_zlevs

     USE mo_param1_bgc,     ONLY : icalc, iopal, idet, idust

     USE mo_memory_bgc,     ONLY : ropal

     USE mo_memory_agg, ONLY     : AJ1, AJ2, AJ3, BJ1, BJ2, BJ3, &
          &                        det_mol2mass, rho_tep, agg_org_dens, &
          &                        dp_dust, dp_det, dp_calc, dp_opal, &
          &                        stickiness_tep, stickiness_det,    &
          &                        stickiness_opal, stickiness_calc,    &
          &                        stickiness_dust, agg_df_min, &
          &                        agg_df_max, agg_re_crit, &
          &                        kavdp, kavrhop, kdfagg, ksticka, klmaxagg, &
          &                        kavrhof
!          &                        kstickf, kbagg, kwsagg, kdynvis, &
!          &                        av_d_c, av_por_V, kavdc, kavpor

     USE mo_hamocc_nml,     ONLY : l_virtual_tep, l_re

     USE mo_sedmnt,         ONLY : claydens, calcdens, opaldens,&
          &                        opalwei, calcwei
 
     USE mo_bgc_memory_types, ONLY  : t_bgc_memory, t_sediment_memory, t_aggregates_memory

     IMPLICIT NONE

     PRIVATE

     PUBLIC:: mean_aggregate_sinking_speed, init_aggregate_params


     REAL(wp), PARAMETER :: one6th = 1._wp/6._wp
     REAL(wp), PARAMETER :: num_fac = 1.e9_wp             ! factor to avoid numerical precision problems
     REAL(wp), PARAMETER :: eps_one = epsilon(1._wp)      ! 2.220446049250313E-016 (Mistral compute2)

     REAL(wp), PARAMETER, DIMENSION(14) :: cdynv=(/1.79e-2_wp, -6.1299e-4_wp,1.4467e-5_wp, &  ! dynamical viscosity
                           &      -1.6826e-7_wp, -1.8266e-7_wp, 9.8972e-12_wp, 2.4727e-5_wp,&
                           &      4.8429e-7_wp, -4.7172e-8_wp,7.5986e-10_wp,1.3817e-8_wp,&
                           &      -2.6363e-10_wp,6.3255e-13_wp,-1.2116e-14_wp/)

! Values given in:
! https://academic.oup.com/plankt/article/27/10/1003/1490556
! following:
! Matthaeus, V. W. (1972) Die Viskositaet des Meerwassers. Beitraege zur Meereskunde, 29, 93-107.
! https://www.io-warnemuende.de/tl_files/forschung/beitraege-zur-meereskunde/1972_29_Beitraege_zur_Meereskunde.pdf
!     REAL(wp),PARAMETER, DIMENSION(14) :: cdynv=(/1.79e-2_wp, -6.1299e-4_wp,1.4467e-5_wp, &
!                           &      -1.6826e-7_wp, -1.7913e-6_wp, 9.5182e-10_wp, 2.4727e-5_wp,&
!                           &      4.8429e-7_wp, -4.7172e-8_wp,7.5986e-10_wp,1.3550e-7_wp,&
!                           &      -2.5853e-9_wp,6.0833e-11_wp,-1.1652e-12_wp/)

     REAL(wp) :: V_dp_dust,V_dp_det,V_dp_calc,V_dp_opal   ! volumes of primary particles (L^3)
     REAL(wp) :: A_dp_dust,A_dp_det,A_dp_calc,A_dp_opal   ! surface areas of primary particles (L^2)
     REAL(wp) :: A_dust,A_det,A_calc,A_opal,A_total       ! total surface area of primary particles per unit volume (L^2/L^3)
     REAL(wp) :: stickiness_min, stickiness_max           ! minimum and maximum stickiness of primary particles
     REAL(wp) :: df_slope                                 ! slope for stickiness to fractal dimension mapping
     REAL(wp) :: rho_V_dp_dust,rho_V_dp_det,rho_V_dp_calc,rho_V_dp_opal ! mass of primary particles (M)
     REAL(wp) :: V_frustule_inner                         ! volume of hollow part in diatom frustule (L^3)
     REAL(wp) :: V_frustule_opal                          ! volume of opal shell material (L^3)
     REAL(wp) :: rho_V_frustule_opal                      ! mass of frustule material (M) 

   CONTAINS

  !=====================================================================================

  SUBROUTINE init_aggregate_params
     !>
     !! Initilization of parameters
     !!
     IMPLICIT NONE
     REAL(wp) :: Rm_SiP                                   ! molar mass ratio opal (SiO_2) to POM
     REAL(wp) :: thick_shell                              ! diatom frustule shell thickness (L)
     REAL(wp) :: d_frustule_inner                         ! diameter of hollow part in diatom frustule (L)
      
     V_dp_dust = one6th * pi * dp_dust**3._wp * num_fac
     V_dp_det  = one6th * pi * dp_det**3._wp  * num_fac
     V_dp_calc = one6th * pi * dp_calc**3._wp * num_fac
     V_dp_opal = one6th * pi * dp_opal**3._wp * num_fac

     A_dp_dust =  pi * dp_dust**2._wp * num_fac
     A_dp_det  =  pi * dp_det**2._wp  * num_fac
     A_dp_calc =  pi * dp_calc**2._wp * num_fac
     A_dp_opal =  pi * dp_opal**2._wp * num_fac

     rho_V_dp_dust = V_dp_dust * claydens
     rho_V_dp_det  = V_dp_det  * agg_org_dens ! orgdens could lead to negative speeds
     rho_V_dp_calc = V_dp_calc * calcdens
     rho_V_dp_opal = V_dp_opal * opaldens

     Rm_SiP          = ropal * opalwei / det_mol2mass
     ! shell thickness
     thick_shell         = 0.5_wp * dp_opal * (1._wp-(opaldens/(Rm_SiP*agg_org_dens+opaldens))**(1._wp/3._wp))
     d_frustule_inner   = dp_opal - 2._wp*thick_shell
     ! volume of hollow part of frustule
     V_frustule_inner    = one6th *  pi  * d_frustule_inner**3._wp * num_fac
     ! volume of opal part of frustule
     V_frustule_opal = one6th *  pi * (dp_opal**3._wp - d_frustule_inner**3._wp) * num_fac
     rho_V_frustule_opal = V_frustule_opal * opaldens
     stickiness_min = MIN(stickiness_tep, stickiness_det,stickiness_opal,stickiness_calc,stickiness_dust)
     stickiness_max = MAX(stickiness_tep, stickiness_det,stickiness_opal,stickiness_calc,stickiness_dust)
     df_slope       = LOG( agg_df_min / agg_df_max)

  END SUBROUTINE init_aggregate_params

  !=====================================================================================

  SUBROUTINE mean_aggregate_sinking_speed (local_bgc_mem, aggr_mem, klev,start_idx, end_idx,pddpo, &
                                       ppo, ptho, psao)
!  SUBROUTINE mean_aggregate_sinking_speed (klev,start_idx, end_idx,pddpo, &
!                                       ppo, ptho, psao, lon, lat)
     TYPE(t_bgc_memory), POINTER    :: local_bgc_mem
     TYPE(t_aggregates_memory), POINTER :: aggr_mem
     
     !-----------------------------------------------------------------------
     !>
     !! calculates the mass concentration-weighted mean sinking velocity of marine
     !! aggregates
     !!

     INTEGER  :: j,k,kpke

     INTEGER, INTENT(in), TARGET    :: klev(bgc_nproma)       !<  vertical levels
     INTEGER, INTENT(in)            :: start_idx              !< start index for j loop (ICON cells, MPIOM lat dir)
     INTEGER, INTENT(in)            :: end_idx                !< end index  for j loop  (ICON cells, MPIOM lat dir)

     REAL(wp), INTENT(in), TARGET   :: pddpo(bgc_nproma,bgc_zlevs)      !< size of scalar grid cell (3rd dimension) [m]
     REAL(wp), INTENT(in), TARGET   :: ppo(bgc_nproma,bgc_zlevs)        !< ocean pressure
     REAL(wp), INTENT(in), TARGET   :: ptho(bgc_nproma,bgc_zlevs)       !< ocean potential temperature [degC]
     REAL(wp), INTENT(in), TARGET   :: psao(bgc_nproma,bgc_zlevs)       !< ocean salinity

!     REAL(wp), INTENT(in)           :: lon(bgc_nproma), lat(bgc_nproma)

     ! molecular dynamic viscosity
     CALL calc_dynvis(aggr_mem%dynvis, klev, start_idx, end_idx, pddpo, ppo, ptho, psao)

     CALL aggregate_properties(local_bgc_mem, aggr_mem, klev, start_idx, end_idx, pddpo, ptho)

     ! calculate the mean sinking velocity of aggregates
     CALL ws_Re_approx(aggr_mem, klev, start_idx, end_idx, pddpo)

     DO j = start_idx, end_idx
        kpke=klev(j)
        DO k = 1,kpke
           IF(pddpo(j,k) > EPSILON(0.5_wp))THEN

!              aggdiag(j,k,kwsagg)    = ws_agg(j,k) * 86400._wp / dtbgc ! conversion  m/time_step   to  m/d for output
!              aggdiag(j,k,kdynvis)   = dynvis(j,k)

              local_bgc_mem%wpoc(j,k)  = aggr_mem%ws_agg(j,k)

              ! set max sinking speed to courant number=0.8, so 0.8*pddpo
              !wpoc(j,k)  = min(wpoc(j,k),0.8_wp*pddpo(j,k))

              local_bgc_mem%wcal(j,k)  = local_bgc_mem%wpoc(j,k)
              local_bgc_mem%wopal(j,k) = local_bgc_mem%wpoc(j,k)
              local_bgc_mem%wdust(j,k) = local_bgc_mem%wpoc(j,k)
           ENDIF
        ENDDO
     ENDDO

  END SUBROUTINE mean_aggregate_sinking_speed

  !=====================================================================================

  SUBROUTINE aggregate_properties(local_bgc_mem, aggr_mem, klev, start_idx, end_idx, pddpo, ptho)
     !-----------------------------------------------------------------------
     !>
     !! aggregate_properties calculates
     !!   - mean stickiness/aggrega
     !!   - fractal dimension
     !!   - slope of aggregate spectrum
     !!   - mean primary particle diameter
     !!   - mean primary particle density
     !!   - maximum aggregate diameter
     !!
     IMPLICIT NONE

     REAL(wp) :: n_det,n_opal,n_calc,n_dust,n_total       ! total primary particle number (#)
     REAL(wp) :: mf                                       ! mass factor for aggregates 
     REAL(wp) :: V_det,V_opal,V_calc,V_dust,V_solid       ! total volume of primary particles in a unit volume (L^3/L^3)
     REAL(wp) :: stickiness_mapped                        ! mapped mean stickiness of particles on range (0,1)
     REAL(wp) :: cell_det_mass                            ! mass of detritus material in diatoms
     REAL(wp) :: cell_pot_det_mass                        ! potential (max) mass detritus material in diatoms
     REAL(wp) :: free_detritus                            ! freely available detritus mass outside the frustule
     REAL(wp) :: V_POM_cell                               ! volume of POM in frustule
     REAL(wp) :: V_aq                                     ! volume of water space in frustule
     REAL(wp) :: rho_frustule                             ! density of diatom frustule incl. opal, detritus and water
     REAL(wp) :: rho_diatom                               ! ensity of either hollow frustule 
 
     TYPE(t_bgc_memory), POINTER    :: local_bgc_mem
     TYPE(t_aggregates_memory), POINTER :: aggr_mem

     INTEGER, INTENT(in), TARGET    :: klev(bgc_nproma)       !<  vertical levels
     INTEGER, INTENT(in)            :: start_idx              !< start index for j loop (ICON cells, MPIOM lat dir)
     INTEGER, INTENT(in)            :: end_idx                !< end index  for j loop  (ICON cells, MPIOM lat dir)

     REAL(wp), INTENT(in), TARGET   :: pddpo(bgc_nproma,bgc_zlevs)      !< size of scalar grid cell (3rd dimension) [m]
     REAL(wp), INTENT(in), TARGET   :: ptho(bgc_nproma,bgc_zlevs)       !< ocean potential temperature [degC]
 
     INTEGER  :: j,k,kpke
     
     REAL(wp),POINTER :: av_dp(:,:),               &  ! mean primary particle diameter
                       &  av_rho_p(:,:),           &  ! mean primary particle density
                       &  df_agg(:,:),             &  ! fractal dimension of aggregates
                       &  b_agg(:,:),              &  ! aggregate number distribution slope
                       &  Lmax_agg(:,:),           &  ! maximum diamater of aggregates
                       &  ws_agg(:,:),             &  ! aggregate mean sinking velocity
                       &  stickiness_agg(:,:),     &  ! mean aggregate stickiness
                       &  stickiness_frustule(:,:),&  ! frustule stickiness
                       &  dynvis(:,:),             &  ! molecular dynamic viscosity
                       &  av_rhof_V(:,:)
     REAL(wp), POINTER :: aggdiag(:,:,:)    ! 3d concentration EU

       av_dp              =>  aggr_mem%av_dp              
       av_rho_p           =>  aggr_mem%av_rho_p           
       df_agg             =>  aggr_mem%df_agg             
       b_agg              =>  aggr_mem%b_agg              
       Lmax_agg           =>  aggr_mem%Lmax_agg           
       ws_agg             =>  aggr_mem%ws_agg             
       stickiness_agg     =>  aggr_mem%stickiness_agg     
       stickiness_frustule=>  aggr_mem%stickiness_frustule
       dynvis             =>  aggr_mem%dynvis             
       av_rhof_V          =>  aggr_mem%av_rhof_V
       aggdiag            =>  aggr_mem%aggdiag

     DO j = start_idx, end_idx
        kpke=klev(j)
        IF(kpke > 0)THEN

        DO k = 1,kpke

         IF(pddpo(j,k) > EPSILON(0.5_wp))THEN

              n_det   = 0._wp ! number of primary particles
              n_opal  = 0._wp
              n_dust  = 0._wp
              n_calc  = 0._wp
              n_total = 0._wp
              mf      = 0._wp

              V_det   = 0._wp ! total volume of primary particles in a unit volume
              V_opal  = 0._wp
              V_calc  = 0._wp
              V_dust  = 0._wp
              V_solid = 0._wp

              free_detritus = 0._wp
              rho_diatom    = 0._wp
              ! n_opal are number of frustule-like primary particles possessing
              ! a density i) different from pure opal ii) due to a mixture of
              ! opal frustule, detritus inside the frustule and potentially water
              ! inside the frustule
              ! that is completely or partially filled with detritus
              ! and water

              cell_det_mass     = 0._wp
              cell_pot_det_mass = 0._wp
              V_POM_cell        = 0._wp
              V_aq              = 0._wp
              rho_frustule      = 0._wp

              ! number of opal frustules (/num_fac)
              n_opal = ABS(local_bgc_mem%bgctra(j,k,iopal)) * opalwei      / rho_V_frustule_opal
! n_opal = 3.124203669937888E-050

              ! maximum mass of detritus inside a frustule
              cell_pot_det_mass = n_opal * V_frustule_inner * agg_org_dens
! cell_pot_det_mass = 1.558762742723002E-052

              ! detritus mass inside frustules
!              cell_det_mass = MIN(cell_pot_det_mass, ABS(bgctra(j,k,idet))  * det_mol2mass - eps_one)
              cell_det_mass = MAX(0._wp,MIN(cell_pot_det_mass, ABS(local_bgc_mem%bgctra(j,k,idet))  * det_mol2mass))
! Examplified values for version "- eps_one"; written out by the write statement at the end of this loop
! cell_det_mass = -2.220446049250313E-016

              ! volume of detritus component in cell
              !              V_POM_cell =   (cell_det_mass/n_opal)/agg_org_dens
              V_POM_cell =   (cell_det_mass/(n_opal+eps_one))/agg_org_dens

              ! if not detritus is available, water is added
              V_aq       =   V_frustule_inner -  V_POM_cell

              ! density of the diatom frsutules incl. opal, detritus and water
              !              rho_frustule = (rho_V_frustule_opal + cell_det_mass/n_opal + V_aq*rhoref_water)/V_dp_opal
              rho_frustule = (rho_V_frustule_opal + cell_det_mass/(n_opal+eps_one) + V_aq*rhoref_water)/V_dp_opal
! rho_frustule = -15027

              ! mass of extra cellular detritus particles
              free_detritus =  ABS(local_bgc_mem%bgctra(j,k,idet))  * det_mol2mass  - cell_det_mass

              IF(l_virtual_tep)then
!                 rho_diatom = (rho_frustule + cell_det_mass/cell_pot_det_mass*rho_tep) &
!                                / (1._wp + cell_det_mass/cell_pot_det_mass)
                 rho_diatom = (rho_frustule + cell_det_mass/(cell_pot_det_mass+eps_one)*rho_tep) &
                                / (1._wp + cell_det_mass/(cell_pot_det_mass+eps_one))
! rho_diatom = -inf; negative cell_det_mass leads to division by 0 here
              ELSE
                 rho_diatom    = rho_frustule
              ENDIF
! rho_diatom = -inf

              ! number of primary particles
              n_det  = free_detritus / rho_V_dp_det  ! includes num_fac
              n_calc = ABS(local_bgc_mem%bgctra(j,k,icalc)) * calcwei      / rho_V_dp_calc
              n_dust = ABS(local_bgc_mem%bgctra(j,k,idust))               / rho_V_dp_dust     ! dust is in kg/m3

               ! total number of primary particles  ! mm: n_total not in use
              ! n_total = n_det + n_opal + n_calc + n_dust

              ! primary particles surface weighted stickiness is mapped
              ! on range between 0 and 1
              ! fractal dimension of aggregates is based on that mapped df
              ! number distribution slope b is based on df

               ! calc total areas
               A_det   = n_det  * A_dp_det
               A_opal  = n_opal * A_dp_opal
               A_calc  = n_calc * A_dp_calc
               A_dust  = n_dust * A_dp_dust
               A_total = A_det + A_opal + A_calc + A_dust

               ! calc frustule stickiness
               stickiness_frustule(j,k) = cell_det_mass / (cell_pot_det_mass + eps_one) * stickiness_tep &
                                          & + (1._wp - cell_det_mass / (cell_pot_det_mass + eps_one)) * stickiness_opal

               ! calc mean stickiness
               stickiness_agg(j,k) = stickiness_frustule(j,k) * A_opal  &
                                     & + stickiness_det  * A_det              &
                                     & + stickiness_calc * A_calc             &
                                     & + stickiness_dust * A_dust

               stickiness_agg(j,k) = stickiness_agg(j,k)/(A_total+eps_one)

               stickiness_mapped     = (stickiness_agg(j,k) - stickiness_min) &
                                     & /(stickiness_max - stickiness_min)

               df_agg(j,k)         = agg_df_max * EXP(df_slope * stickiness_mapped)

               ! Slope is here positive defined (as n(d)~d^-b), so *-1 of
               ! Jiang & Logan 1991: Fractal dimensions of aggregates
               ! determined from steady-state size distributions.
               ! Environ. Sci. Technol. 25, 2031-2038.
               !
               ! See also:
               ! Hunt 1980: Prediction of oceanic particle size distributions
               !            from coagulation and sedimentation mechanisms.
               !
               ! Additional assumptions made here:
               ! b in Jiang & Logan     (used for       Re <   0.1: b=1
               !                              for 0.1 < Re <  10  : b=0.871
               !                              for 10  < Re < 100  : b=0.547)
               ! is set to 0.871 as an 'average for our range of 0<Re<Re_crit'
               ! D2=min(2,df(3d)) (Meakin 1988)
               !
               ! => Formulation in Jiang & Logan 1991:
               ! slope = -0.5*(3+df+(2+df-D2)/(2-b)) reduces to:

               b_agg(j,k) = 0.5_wp * (3._wp + df_agg(j,k) &
                         & + (2._wp + df_agg(j,k) - &
                         & MIN(2._wp, df_agg(j,k)) ) / (2._wp - BJ2))

               ! careful: for df=1.5904: b_agg=2*df where w_s is undefined.

               ! total volume of primary particles
               V_det   = n_det  * V_dp_det       * num_fac
               V_opal  = n_opal * V_dp_opal      * num_fac
               V_calc  = n_calc * V_dp_calc      * num_fac
               V_dust  = n_dust * V_dp_dust      * num_fac
               V_solid = V_det + V_opal + V_calc + V_dust

               ! primary particle mean diameter according to Bushell & Amal 1998, 2000
               ! sum(n_i) not changing - can be pulled out and thus cancels out
               av_dp(j,k) = &
                 & (n_calc*dp_calc**3._wp + n_dust*dp_dust**3._wp + n_opal*dp_opal**3._wp + n_det*dp_det**3._wp)
               
               av_dp(j,k) = av_dp(j,k) / &
                 & (n_calc*dp_calc**df_agg(j,k) + n_dust*dp_dust**df_agg(j,k) &
                 & + n_opal*dp_opal**df_agg(j,k) + n_det*dp_det**df_agg(j,k) + eps_one)
               
               av_dp(j,k) = av_dp(j,k)**(1._wp/(3._wp-df_agg(j,k)))

               ! density of mean primary particles
               av_rho_p(j,k) = (V_det*agg_org_dens + V_opal*rho_diatom + V_calc*calcdens + V_dust*claydens) &
                                  & / (V_solid + eps_one)
!               if (v_solid == 0._wp) write(*,*) '#1# j, k, v_solid, av_rho_p:',j,k,v_solid,av_rho_p(j,k)
!               if (av_rho_p(j,k) < 0._wp .or. av_rho_p(j,k) > 2600._wp) then
!                  write(*,*) '#1#:',j,k,av_rho_p(j,k),v_det,v_opal,v_calc,v_dust,v_solid,&
!                       & rho_diatom, rho_frustule, cell_det_mass, cell_pot_det_mass, rho_tep,&
!                       & n_opal, V_frustule_inner, agg_org_dens,&
!                       & bgctra(j,k,idet), det_mol2mass, eps_one
!               endif

            ENDIF
         ENDDO
         ENDIF
      ENDDO
 
 

   ! calculate the maximum diameter of aggregates based on agg props
!   CALL max_agg_diam(kpie, kpje, kpke, pddpo)
   CALL max_agg_diam(aggr_mem, klev, start_idx, end_idx, pddpo)

 
     DO j = start_idx, end_idx
        kpke=klev(j)
        IF(kpke > 0)THEN

        DO k = 1,kpke

         IF(pddpo(j,k) > EPSILON(0.5_wp))THEN

             ! mass factor  ! mm: only used to calculate number of aggregates n_agg; see unmodularized MAGO version
!             mf = mass_factor(av_dp(j,k), df_agg(j,k), av_rho_p(j,k))

!             av_d_c(j,k)    = (1._wp + df_agg(j,k) - b_agg(j,k))                    &
!                            & /(2._wp + df_agg(j,k) - b_agg(j,k))                   &
!                            & *(Lmax_agg(j,k)**(2._wp + df_agg(j,k) - b_agg(j,k))   &
!                            & - av_dp(j,k)**(2._wp + df_agg(j,k) - b_agg(j,k)))     &
!                            & / (Lmax_agg(j,k)**(1._wp+df_agg(j,k)-b_agg(j,k))      &
!                            & - av_dp(j,k)**(1._wp + df_agg(j,k)-b_agg(j,k)))
!
             ! volume-weighted aggregate density
             av_rhof_V(j,k) = &
              & (av_rho_p(j,k) - rhoref_water) &
              & * av_dp(j,k)**(3._wp - df_agg(j,k)) &
              & * (4._wp - b_agg(j,k)) * &
              & (Lmax_agg(j,k)**(1._wp + df_agg(j,k) - b_agg(j,k)) &
              & - av_dp(j,k)**(1._wp + df_agg(j,k) - b_agg(j,k)))  &
              & / ((1._wp + df_agg(j,k) - b_agg(j,k))                       &
              & * (Lmax_agg(j,k)**(4._wp - b_agg(j,k)) - &
              &    av_dp(j,k)**(4._wp - b_agg(j,k))))   &
              & + rhoref_water
!
!             ! volume-weighted aggregate porosity
!             av_por_V(j,k)  =  1._wp - ((4._wp - b_agg(j,k))                       &
!                            & * av_dp(j,k)**(3._wp - df_agg(j,k))                  &
!                            & * (Lmax_agg(j,k)**(1._wp + df_agg(j,k) - b_agg(j,k)) &
!                            & - av_dp(j,k)**(1._wp + df_agg(j,k) - b_agg(j,k))))   &
!                            & / ((1._wp + df_agg(j,k) - b_agg(j,k))                &
!                            & * (Lmax_agg(j,k)**(4._wp - b_agg(j,k)) - av_dp(j,k)**(4._wp - b_agg(j,k))))

             aggdiag(j,k,kdfagg)    = df_agg(j,k)
             aggdiag(j,k,klmaxagg)  = lmax_agg(j,k)
             aggdiag(j,k,kavdp)     = av_dp(j,k)
             aggdiag(j,k,ksticka)    = stickiness_agg(j,k)
             aggdiag(j,k,kavrhop)   = av_rho_p(j,k)
             aggdiag(j,k,kavrhof)   = av_rhof_V(j,k)
!             aggdiag(j,k,kstickf)   = stickiness_frustule(j,k)
!             aggdiag(j,k,kavdc)     = av_d_c(j,k)
!             aggdiag(j,k,kbagg)     = b_agg(j,k)
!             aggdiag(j,k,kavpor)    = av_por_V(j,k)

            ENDIF
         ENDDO
         ENDIF
      ENDDO
 
 

  END SUBROUTINE aggregate_properties


  !=====================================================================================

  SUBROUTINE ws_Re_approx(aggr_mem, klev, start_idx, end_idx, pddpo)
     !-----------------------------------------------------------------------
     !>
     !! ws_Re_approx:  distribution integrated to Lmax (Re crit dependent maximum agg size)
     !! Renolds number-dependent sinking velocity.
     !! Approximation for c_D-value taken from Jiang & Logan 1991:
     !! c_D=a*Re^-b
     !!
     IMPLICIT NONE

     TYPE(t_aggregates_memory), POINTER :: aggr_mem
     INTEGER  :: j,k,kpke

     INTEGER, INTENT(in), TARGET    :: klev(bgc_nproma)       !<  vertical levels
     INTEGER, INTENT(in)            :: start_idx              !< start index for j loop (ICON cells, MPIOM lat dir)
     INTEGER, INTENT(in)            :: end_idx                !< end index  for j loop  (ICON cells, MPIOM lat dir)
     REAL(wp), INTENT(in), TARGET   :: pddpo(bgc_nproma,bgc_zlevs)      !< size of scalar grid cell (3rd dimension) [m]

 
     DO j = start_idx, end_idx
        kpke=klev(j)
        IF(kpke > 0)THEN
        DO k = 1,kpke
           IF(pddpo(j,k) > EPSILON(0.5_wp))THEN
              ! ws_Re is a function
              aggr_mem%ws_agg(j,k) = ws_Re(aggr_mem,j,k)
           ENDIF
        ENDDO
        ENDIF
     ENDDO
 
 

  END SUBROUTINE ws_Re_approx

  !=====================================================================================

  REAL(wp) FUNCTION get_dRe(dynvis, df_agg, av_rho_p, av_dp,  AJ, BJ, Re)
     IMPLICIT NONE
     ! Arguments
     REAL(wp) :: dynvis, df_agg, av_rho_p, av_dp
     REAL(wp), INTENT(in) :: AJ
     REAL(wp), INTENT(in) :: BJ
     REAL(wp), INTENT(in) :: Re

    ! Local variables
    REAL(wp) :: nu_vis

    nu_vis =  dynvis/rhoref_water
    get_dRe = (Re*nu_vis)**((2._wp - BJ)/df_agg)/(4._wp/3._wp*(av_rho_p - rhoref_water)/rhoref_water &
           *av_dp**(3._wp - df_agg)*g/(AJ*nu_vis**(BJ)))**(1._wp/df_agg)

  END FUNCTION get_dRe

  REAL(wp) FUNCTION get_ws_agg_integral(dynvis, av_rho_p, av_dp, df_agg, b_agg, AJ, BJ, lower_bound, upper_bound)
    IMPLICIT NONE

    REAL(wp) :: dynvis, av_rho_p, av_dp, df_agg, b_agg
    REAL(wp), INTENT(in) :: AJ
    REAL(wp), INTENT(in) :: BJ
    REAL(wp), INTENT(in) :: upper_bound
    REAL(wp), INTENT(in) :: lower_bound

    ! Local variables
    REAL(wp) :: nu_vis

    nu_vis = dynvis/rhoref_water
    get_ws_agg_integral = (4._wp/3._wp*(av_rho_p - rhoref_water)/rhoref_water &
                     & *av_dp**(3._wp - df_agg)*g  &
                     & /(AJ*nu_vis**BJ))**(1._wp/(2._wp - BJ)) &
                     & *(upper_bound**(1._wp - b_agg + df_agg    &
                     & + (BJ + df_agg - 2._wp)/(2._wp - BJ)) &
                     & /(1._wp - b_agg + df_agg + (BJ + df_agg - 2._wp)/(2._wp - BJ)) &
                     & - lower_bound**(1._wp - b_agg + df_agg + (BJ + df_agg - 2._wp)   &
                     & /(2._wp - BJ)) &
                     & /(1._wp - b_agg + df_agg + (BJ + df_agg - 2._wp)/(2._wp - BJ)))

  END FUNCTION get_ws_agg_integral

  REAL(wp) FUNCTION ws_Re(aggr_mem, j,k)
    !-----------------------------------------------------------------------
    !>
    !! ws_Re:  distribution integrated to Lmax (Re crit dependent maximum agg size)
    !! Reynolds number-dependent sinking velocity.
    !! Approximation for c_D-value taken from Jiang & Logan 1991:
    !! c_D=a*Re^-b
    !! written in such a way that we check the critical Reynolds
    !! number (in case that we extend the maximum size by shear-
    !! driven break-up).
    !!

    IMPLICIT NONE

    TYPE(t_aggregates_memory), POINTER :: aggr_mem
    INTEGER, INTENT(in)  :: j
    INTEGER, INTENT(in)  :: k

    ! Local
    REAL(wp) :: d_Re01, d_Re10, d_low, ws_agg_ints
    REAL(wp),POINTER :: av_dp(:,:),               &  ! mean primary particle diameter
                       &  av_rho_p(:,:),           &  ! mean primary particle density
                       &  df_agg(:,:),             &  ! fractal dimension of aggregates
                       &  b_agg(:,:),              &  ! aggregate number distribution slope
                       &  Lmax_agg(:,:),           &  ! maximum diamater of aggregates
                       &  ws_agg(:,:),             &  ! aggregate mean sinking velocity
                       &  stickiness_agg(:,:),     &  ! mean aggregate stickiness
                       &  stickiness_frustule(:,:),&  ! frustule stickiness
                       &  dynvis(:,:),             &  ! molecular dynamic viscosity
                       &  av_rhof_V(:,:)

       av_dp              =>  aggr_mem%av_dp              
       av_rho_p           =>  aggr_mem%av_rho_p           
       df_agg             =>  aggr_mem%df_agg             
       b_agg              =>  aggr_mem%b_agg              
       Lmax_agg           =>  aggr_mem%Lmax_agg           
       ws_agg             =>  aggr_mem%ws_agg             
       stickiness_agg     =>  aggr_mem%stickiness_agg     
       stickiness_frustule=>  aggr_mem%stickiness_frustule
       dynvis             =>  aggr_mem%dynvis             
       av_rhof_V          =>  aggr_mem%av_rhof_V

    ! for Re-dependent, it should always be agg_Re_crit>10
    ! for shear-driven break-up, check against integration bounds
    ! calc integration limits for Re-dependent sinking:
    ! Re=0.1
    d_Re01 = get_dRe(dynvis(j,k), df_agg(j,k), av_rho_p(j,k), av_dp(j,k), &
      & AJ1, BJ1, 0.1_wp)
    
    ! Re=10
    d_Re10 = get_dRe(dynvis(j,k), df_agg(j,k), av_rho_p(j,k), av_dp(j,k), &
      & AJ2, BJ2, 10._wp)
    d_low = av_dp(j,k)

    ws_agg_ints = 0._wp
    IF(Lmax_agg(j,k) >= d_Re01)THEN ! Re > 0.1
                                    ! - collect full range up to
                                    ! 0.1, (dp->d_Re1) and set lower bound to
                                    ! Re=0.1 val
                                    ! aj=AJ1, bj=1
      ws_agg_ints = get_ws_agg_integral(dynvis(j,k), av_rho_p(j,k), av_dp(j,k), df_agg(j,k), b_agg(j,k), &
        & AJ1, BJ1, av_dp(j,k), d_Re01)
        d_low = d_Re01
    ENDIF

    IF(Lmax_agg(j,k) >= d_Re10)THEN ! Re > 10
                                         ! - collect full range Re=0.1-10 (d_Re1-> d_Re2)
                                         ! and set lower bound to
                                         ! Re=10 val
                                         ! aj=AJ2, bj=0.871
        ws_agg_ints = ws_agg_ints  + &
         get_ws_agg_integral(dynvis(j,k), av_rho_p(j,k), av_dp(j,k), df_agg(j,k), b_agg(j,k), &
          AJ2, BJ2, d_Re01, d_Re10)
        d_low = d_Re10
    ENDIF

    IF(d_low < d_Re01)THEN ! Re<0.1 and Lmax < d_Re1
        ws_agg_ints = get_ws_agg_integral(dynvis(j,k), av_rho_p(j,k), av_dp(j,k), df_agg(j,k), b_agg(j,k), &
          & AJ1, BJ1, av_dp(j,k), Lmax_agg(j,k))
    ELSE ! Re > 10, aj=AJ3, bj=BJ3
        ws_agg_ints = ws_agg_ints + &
          & get_ws_agg_integral(dynvis(j,k), av_rho_p(j,k), av_dp(j,k), df_agg(j,k), b_agg(j,k), AJ3, BJ3, d_low, Lmax_agg(j,k))
    ENDIF

    ! concentration-weighted mean sinking velocity
    ws_Re = (ws_agg_ints &
            & /((Lmax_agg(j,k)**(1._wp + df_agg(j,k) - b_agg(j,k))  &
            & - av_dp(j,k)**(1._wp + df_agg(j,k) - b_agg(j,k)))  &
            & / (1._wp + df_agg(j,k) - b_agg(j,k))))*dtbgc   ! (m/s -> m/d)  *dtb

  END FUNCTION ws_Re

  !=====================================================================================

  SUBROUTINE max_agg_diam(aggr_mem, klev, start_idx, end_idx, pddpo)
     !-----------------------------------------------------------------------
     !>
     !! max_agg_diam calculates the maximum aggregate diameter of the aggregate
     !! number distribution, assumes Re_crit > 10
     !!

     TYPE(t_aggregates_memory), POINTER :: aggr_mem
     
     INTEGER  :: j,k,kpke

     INTEGER, INTENT(in), TARGET    :: klev(bgc_nproma)       !<  vertical levels
     INTEGER, INTENT(in)            :: start_idx              !< start index for j loop (ICON cells, MPIOM lat dir)
     INTEGER, INTENT(in)            :: end_idx                !< end index  for j loop  (ICON cells, MPIOM lat dir)

     REAL(wp), INTENT(in), TARGET   :: pddpo(bgc_nproma,bgc_zlevs)      !< size of scalar grid cell (3rd dimension) [m]

     ! Local variables
     REAL(wp) :: nu_vis

 
     DO j = start_idx, end_idx
        kpke=klev(j)
        IF(kpke > 0)THEN

        DO k = 1,kpke

           IF(pddpo(j,k) > EPSILON(0.5_wp))THEN

! ori:
!              Lmax_agg(j,k) = ((agg_Re_crit * 18._wp *  dynvis(j,k) * dynvis(j,k) / rhoref_water)&
!                                 & /((av_rho_p - rhoref_water)*g))**(1._wp/df_agg) &
!                                 & * av_dp**(1._wp-3._wp/df_agg)

! new: based on analytical Jiang approximation
              nu_vis  =  aggr_mem%dynvis(j,k)/rhoref_water
              aggr_mem%Lmax_agg(j,k) = (agg_Re_crit*nu_vis)**((2._wp - BJ3)/aggr_mem%df_agg(j,k))        &
                                & /((4._wp/3._wp)*(aggr_mem%av_rho_p(j,k) - rhoref_water)/rhoref_water    &
                                & *aggr_mem%av_dp(j,k)**(3._wp - aggr_mem%df_agg(j,k))*g &
                                & /(AJ3*nu_vis**BJ3))**(1._wp/aggr_mem%df_agg(j,k))

           ENDIF
        ENDDO
        ENDIF
     ENDDO
 
 

  END SUBROUTINE max_agg_diam

  !=====================================================================================

  REAL (wp) FUNCTION mass_factor(dp,df,rhop)
     !-----------------------------------------------------------------------
     !>
     !! mass_factor calculates the mass factor for the mass of a single
     !! aggregate
     !!
    IMPLICIT NONE

    REAL(wp), INTENT(in) :: dp
    REAL(wp), INTENT(in) :: df
    REAL(wp), INTENT(in) :: rhop

     ! mass factor
     mass_factor = one6th * pi * dp**(3._wp - df) * rhop

  END FUNCTION mass_factor

  REAL(wp) FUNCTION rho_agg(d,rhop,dp,df,rho)
     !-----------------------------------------------------------------------
     !>
     !! rho_agg provides the aggregate density
     !!

     IMPLICIT NONE

     REAL(wp), INTENT(in) :: d
     REAL(wp), INTENT(in) :: rhop
     REAL(wp), INTENT(in) :: dp
     REAL(wp), INTENT(in) :: df
     REAL(wp), INTENT(in) :: rho

     rho_agg =  (rhop - rho)*(dp/d)**(3._wp - df) + rho

  END FUNCTION rho_agg

  !=====================================================================================

  SUBROUTINE calc_dynvis(dynvis, klev, start_idx, end_idx, pddpo, ppo, ptho, psao)
     !-----------------------------------------------------------------------
     !>
     !! dynvis calculates the molecular dynamic viscosity according to
     !! Richards 1998: The effect of temperature, pressure, and salinity
     !! on sound attenuation in turbid seawater. J. Acoust. Soc. Am. 103 (1),
     !! originally published by  Matthaeus, W. (1972): Die Viskositaet des
     !! Meerwassers. Beitraege zur Meereskunde, Heft 29 (in German).
     !!

     IMPLICIT NONE
     REAL(wp) :: dynvis(:,:)
     
     INTEGER  :: j,k,kpke

     INTEGER, INTENT(in), TARGET    :: klev(bgc_nproma)       !<  vertical levels
     INTEGER, INTENT(in)            :: start_idx              !< start index for j loop (ICON cells, MPIOM lat dir)
     INTEGER, INTENT(in)            :: end_idx                !< end index  for j loop  (ICON cells, MPIOM lat dir)

     REAL(wp), INTENT(in), TARGET   :: pddpo(bgc_nproma,bgc_zlevs)      !< size of scalar grid cell (3rd dimension) [m]
     REAL(wp), INTENT(in), TARGET   :: ppo(bgc_nproma,bgc_zlevs)        !< ocean pressure
     REAL(wp), INTENT(in), TARGET   :: ptho(bgc_nproma,bgc_zlevs)       !< ocean potential temperature [degC]
     REAL(wp), INTENT(in), TARGET   :: psao(bgc_nproma,bgc_zlevs)       !< ocean salinity

     ! Local variables
     REAL(wp):: press_val  ! Pascal/rho -> dbar 

 
     DO j = start_idx, end_idx
        kpke=klev(j)
        IF(kpke > 0)THEN

        DO k = 1,kpke

           IF(pddpo(j,k) > EPSILON(0.5_wp))THEN

!              press_val = ppo(j,k) * rhoref_water * 1.e-5_wp ! Pascal/rho -> dbar
              press_val = ppo(j,k) * rhoref_water/g * 1.e-4_wp ! mm
! MAGO needs dynvis to be in [Pa*s = kg m-1 s-1], as provided by the formula here.
! The formula needs press_val to be in [kg cm-2]. p=rho*g*h -> p/g is [kg m-2] -> p/g * 1.e-4 is [kg cm-2]
! ppo is provided as p/rho [m2/s2] (see ocean/physics/mo_ocean_thermodyn.f90; p is NOT in [m] as mentioned for press_hyd in ocean/dynamics/mo_ocean_types.f90)
! To get press_val in [kg cm-2], ppo has to be multiplied by rho/g and 1.e-4.

              ! molecular dynamic viscosity
              dynvis(j,k) = 0.1_wp *       & ! Unit: g / (cm*s) -> kg / (m*s)
                   &       (cdynv(1)                                                                           &
                   &      + cdynv(2)*ptho(j,k) + cdynv(3)*ptho(j,k)**2._wp + cdynv(4)*ptho(j,k)**3._wp &
                   &      + cdynv(5)*press_val + cdynv(6)*press_val**2._wp &
                   &      + cdynv(7)*psao(j,k) &
                   &      + psao(j,k)*(cdynv(8)*ptho(j,k) + cdynv(9)*ptho(j,k)**2._wp + cdynv(10)*ptho(j,k)**3._wp) &
                   &      + press_val*(cdynv(11)*ptho(j,k) + cdynv(12)*ptho(j,k)**2._wp)         &
                   &      - press_val**2._wp*(cdynv(13)*ptho(j,k) + cdynv(14)*ptho(j,k)**2._wp) &
                   &        )

           ENDIF
        ENDDO
        ENDIF
     ENDDO
 
 

  END SUBROUTINE calc_dynvis

END MODULE mo_aggregates

