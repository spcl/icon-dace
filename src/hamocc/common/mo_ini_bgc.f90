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

! @brief set start values for bgc variables

MODULE mo_ini_bgc

  USE mo_kind, ONLY        : wp
  USE mo_bgc_memory_types, ONLY  :   t_bgc_memory, t_sediment_memory, t_aggregates_memory



  USE mo_memory_bgc, ONLY   :atmacon, atmacmol,    &
       &                     ozkoa,ralk, ro2ut_cya, cyamin,    &
       &                     calcinp,orginp,silinp, &
       &                     phytomi, grami, remido, dyphy, zinges,        &
       &                     epsher,  spemor, gammap, gammaz, ecan, &
       &                     pi_alpha, fpar, bkphy, bkzoo, bkopal,         &
       &                     dremn2o, sulfate_reduction,         &
       &                     n2_fixation, ro2ut, rcar, rnit,     &
       &                     rnoi, nitdem, n2prod, ropal,   &
       &                     perc_diron, riron, fesoly, relaxfe,     &
       &                     rn2,             &
       &                     thresh_o2,   &
       &                     pi_alpha_cya,          &
       &                     Topt_cya,T1_cya,T2_cya,bkcya_N, &
       &                     buoyancyspeed_cya, bkh2sox, rh2sox, &
       &                     doccya_fac, thresh_aerob, thresh_sred, &
       &                     wcya, p2gtc, ro2bal, dmsp,prodn2o,docmin, &
       &                     no2denit, anamoxra, nitriox, nitrira, ro2ammo, &
       &                     bknh4_cya, bkno3_cya, bkno3, bknh4, rmm, kg_denom, bkpo4, &
       &                     bkno2, bkrad, bkfe, rno3nh4, rno3no2, rno2no3, rnh4no2, &
       &                     alk_nrn2, rno2n2, o2thresh, o2den_lim, rrrcl,  &
       &                     sinkspeed_dust

  USE mo_memory_agg, ONLY  : agg_org_dens, det_mol2mass, rho_tep, &
       &                     AJ1, AJ2, AJ3, BJ1, BJ2, BJ3, &
       &                     dp_dust, dp_det, dp_calc, dp_opal, &
       &                     stickiness_tep, stickiness_det, stickiness_opal, &
       &                     stickiness_calc, stickiness_dust, &
       &                     agg_df_min, agg_df_max, agg_re_crit

  USE mo_sedmnt, ONLY      : disso_op,disso_cal,&
       &                     o2ut, rno3, sred_sed, silsat, calcon

  USE mo_hamocc_nml, ONLY  : l_cpl_co2, i_settling, &
       &                     sinkspeed_poc, sinkspeed_opal, sinkspeed_calc,&
       &                     ks,cycdec,cya_growth_max,grazra,&
       &                     mc_fac, sinkspeed_martin_ez, mc_depth, denit_sed, disso_po, &
       &                     atm_co2, atm_o2, atm_n2, deltacalc, deltaorg, deltasil, &
       &                     drempoc, dremopal, dremcalc,  denitrification, &
       &                     l_N_cycle, no3nh4red, no3no2red


  USE mo_control_bgc, ONLY : ldtbgc, dtb, dtbgc, rmasko, rmasks, &
       &                     bgc_gin,  & 
       &                     bgc_tatl, &
       &                     bgc_land,  &
       &                     bgc_soce, bgc_npac,  &
       &                     bgc_nproma, bgc_zlevs

  USE mo_param1_bgc, ONLY  : ialkali, ian2o, iatmco2, iatmn2,     &
       &                     iano3, iatmo2, icalc, idet, idoc, igasnit,   &
       &                     iopal, ioxygen, idust,iagesc,                &
       &                     iphosph, iphy, ipowaal, ipowaic, ipowaox,    &
       &                     ipowaph, ipowasi, ipown2, ipowno3,           &
       &                     isco212, isilica, isssc12, issso12, issssil, &
       &                     izoo, ipowafe, issster, &
       &                     icya, iiron, idms, ih2s, ipowh2s, &
       &                     iammo, iano2, ipownh4, ipowno2

  USE mo_bgc_constants
     
!  USE mo_planetary_constants, ONLY: g, rhoref_water
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ini_aquatic_tracers,            &
       &    ini_pore_water_tracers,         &
       &    ini_atmospheric_concentrations, &
       &    ini_wpoc, bgc_param_conv_unit,  &
       &    ini_continental_carbon_input,   &
       &    set_parameters_bgc,             & 
       &    ini_aggregate_parameters !,             &
!       &    level_ini





CONTAINS


  SUBROUTINE SET_PARAMETERS_BGC
    !
    ! Initialize overall time step counter.
    !
    ldtbgc = 0

    ozkoa    = 3.e15_wp                      ! ozean kohlenstoff alt, start with dummy value in maschk
    atmacon  = 0.35e-3_wp * 5.1e14_wp*12._wp ! factor for atmospheric concentration of CO2 in g C -> ppm
    atmacmol = 0.35e-3_wp * 5.1e14_wp        ! factor for atmospheric concentration of CO2 in mol C -> ppm
    ! 5.1e14 = total surface of Earth in [m^2]

    !
    ! Biology
    !
    phytomi  = 1.e-11_wp   ! i.e. 1e-5 mmol P/m3 minimum concentration of phytoplankton
    grami    = 1.e-11_wp   ! i.e. 1e-5 mmol P/m3 minimum concentration of zooplankton | test e-11 ->e-8 very slow decay
    cyamin   = 1.e-11_wp   ! minimum cyanobacteria conc. (< initial value)
    docmin   = 1.e-10_wp

    remido   = 0.008_wp     ! KS, JS, EMR 12/07/2007
    dyphy    = 0.008_wp     ! 1/d -mortality rate of phytoplankton
    zinges   = 0.6_wp            ! dimensionless fraction -assimilation efficiency of zooplankton
    epsher   = 0.8_wp            ! dimensionless fraction - (1-epsher)=fraction of grazing egested
    spemor   = 3.e6_wp      ! 1/d -mortality rate of zooplankton
    gammap   = 0.03_wp      ! 1/d -exudation rate
    gammaz   = 0.06_wp      ! 1/d -excretion rate
    ecan     = 0.95_wp           ! fraction of mortality as PO_4
    pi_alpha = 0.02_wp           ! initial slope of production vs irradiance
    fPAR     = 0.4_wp            ! fraction of Photosynthetic Active Radiation
    thresh_aerob = 5.e-8_wp      ! kmol m-3,  O2 threshold for aerob remineralization
    thresh_o2 = 10.e-6_wp      ! kmol m-3,  O2 threshold for aerob remineralization
    thresh_sred = 3.e-6_wp      ! kmol m-3,  O2 threshold for sulfate reduction
    prodn2o = 1.e-4_wp



    bkphy  = 1.e-8_wp 
    bkzoo  = 4.e-8_wp
    bkopal = 1.e-6_wp 
    bkh2sox = 5.e-7_wp ! i.e. 0.5 mmol O2 m-3 for H2S oxidation



    ! water column remineralisation constants

    rh2sox  = 0.93_wp ! 1/d  H2S oxidation rate

! ------ cyanobacteria
    buoyancyspeed_cya = 1._wp   ! daily buoyancy speed of cya  
    pi_alpha_cya      = 0.03_wp      ! m2 W-1 d-1
    Topt_cya          = 28._wp       ! deg C
    T1_cya            = 5.5_wp       ! deg C
    T2_cya            = 1._wp        ! deg C
    bkcya_N           = 1.e-9_wp     ! kmol/m3  
    doccya_fac        = 0.1_wp
! ------

    dremn2o  = 0.01_wp      ! 1/d
    sulfate_reduction = 0.005_wp

    ! nitrogen fixation 
    n2_fixation = 0.005_wp

    ! extended redfield ratio declaration
    ! Note: stoichiometric ratios are based on Takahashi etal. (1985)
    ! P:N:C:-O2 + 1:16:122:172

    ro2ut = 172._wp
    ro2ut_cya = 148._wp ! ro2ut - 3 * rnit/rn2 = 172 - 3*16/2 = 148
    rcar  = 122._wp
    rnit  = 16._wp
    rnoi  = 1._wp/rnit
    rn2   = 2._wp
    ro2bal = ro2ut-rcar-2._wp-1.5_wp*rnit ! OM conversion for O2 mass balance
                                          ! OM { rcar C: rnit N : P }
                                          ! C -> CO2, N-> NO3, P -> PO4
                                          ! C[O2] = 1, N[O2] =1.5, P[O2]=2 
                                          ! total O2   = [O2]+[CO2]+2[PO4]+1.5[NO3] - ro2bal[OM]   
    ralk = rnit + 1._wp                   ! for alkalinity updates during OM prod/loss
                                          ! 16 H+ for NO3, 1 for P assim/release 
    p2gtc = rcar*12._wp*1.e-12_wp ! kmolP to GtC
    

    ! N consumption of denitrification corrected after Paulmier etal, 2009)
    nitdem = 137.6_wp          ! NO3 demand of denitrification
    n2prod = nitdem * 0.5_wp   ! N2 production during denitrification

    ropal  = 25._wp ! opal to organic phosphorous production ratio

    ! for interaction with sediment module
    o2ut = 172._wp
    rno3 = 16._wp

    ! Sediment
    denit_sed = 0.01_wp      ! 1/d  sediment denitrification rate 
    sred_sed  = 0.0001_wp    ! sediment sulfate reduction rate
    disso_op  = 0.002592_wp  ! [m3/(kmol Si(OH)4)*1/d ]dissolution rate of opal 
    disso_cal = 0.00864_wp   ! [m3/(kmol CO3) 1/d] dissolution rate of CaCO3
    disso_po  = 0.01_wp      ! [m3/(kmol O2) 1/d] degradation rate organic C in sediment
    silsat    = 0.001_wp     ! [mol/m3] Silicate saturation constant  
    o2thresh  = 10.E-6_wp    ! set to 10 umol/l O2 threshold in sediment
    o2den_lim = 0.5E-6_wp


    ! weight percent iron in dust deposition (0.035) times Fe solubility (0.01) /55.85 g--> mol
    perc_diron = 0.035_wp * 0.01_wp / 55.85_wp

    ! the three values below are from Johnson et al., 1997 Mar. Chemistry 57, 137-161
    ! riron   = 5.*rcar*1.e-6       ! 'Q'=5.*122.*1.e-6 = 6.1e-4   (iron to phosphate stoichiometric ratio * umol->mol)
    ! riron   = 2.*rcar*1.e-6       ! 14.2.06: 5.->2. (2.5umol P/nmol Fe) emr: 3.*rcar*1.e-6
    riron   = 3._wp * rcar*1.e-6_wp ! 06.3.06: 2.->3. coex90=12.2GTC/a before change

    fesoly  = 0.6_wp *1.e-9_wp      ! global mean/max (?) dissolved iron concentration in deep water (for relaxation) [mol/m3]
    ! 'fesoly' stands for 'apparent solubility value' of iron

    relaxfe = 0.05_wp/365._wp       ! relaxation time for iron to fesoly corresponds to 20 yrs
    !  relaxfe = 0.005_wp/365._wp*dtb       ! changed as suggested by emr
    ! (relaxfe is only used to remove iron, corresponding to Johnson's 200yr residence time)


    ! DMS
     dmsp(1) = 10._wp     ! temperature dependent release by phytoplankton
     dmsp(2) = 0.0075_wp  ! photolysis (uv-decay)
     dmsp(3) = 0.0096_wp  ! microbial consumption
     dmsp(4) = 1.25_wp * 0.107638502e0_wp  ! production with delcar(coccolithoporides) 
     dmsp(5) = 1.25_wp * 0.109784522e-1_wp  ! production with delsil(diatoms)
     dmsp(6) = 0.1e-07_wp            ! half saturation rate const. bacterial decomp
     !$ACC ENTER DATA COPYIN(dmsp)

    !
    ! extended N-cycle
    !
    bkno3 = 0.16_wp*1.e-6_wp
    bknh4 = 0.1_wp *1.e-6_wp
    ro2ammo = ro2ut - 2._wp*rnit
   !       oxygen demand to nitrify nh4 to no2
    rnh4no2 =  24._wp /rnit 
    alk_nrn2 = (560._wp + 48._wp)/3._wp ! alkalinity is increased acc
!      oxygen  demand to nitrify no2 to no3
    rno2no3 = 8._wp /rnit ! oxygen demand during nitrification in P-units 2*rnit, 
    bkno3_cya = bkno3*1.e-1_wp
    bknh4_cya = bknh4*1e-1_wp
    bkpo4 = 0.01_wp*1.E-6_wp   ! in kmolP/m3 half satur. const. for PO4 
    bkfe = bkpo4*riron
    rno2n2 = 560._wp/3._wp 
    rmm = 17.03_wp
    kg_denom = 770._wp + 45._wp*rmm**(1._wp/3._wp)  ! denominator of gas-phase tranfser velocity
!                                               for ammonia with 17.03 relative mol. mass of NH3
!   Chemical ratios for constituents  during remineralization,
!   denitrification and so on (after beckmann&hense ,2012, adapted 
! to our Redfield ratios

    ! DNRN : NO3 reduction to NO2
      rno3no2  = 280._wp  ! nitrate used and  nitrite produced per P-unit org

    ! DNRA : NO3 reduction to NH4
      rno3nh4  = 70._wp   ! nitrate used per P-unit org

    ! NRN2   :  NO2 to N2
      no2denit = 0.008_wp ! 1/day 
   
    !ANAMMOX
      anamoxra = 0.05_wp         ! anammox rate  1/day 
      bkno2 = 0.5_wp*1.E-6_wp    ! Half saturation constant for Nitrite in kmolN/m3
    !NITOX : oxidation of NO2 to NO3 ; light dependent	   
      nitriox = 0.25_wp          !   nitrite oxidation rate 1/day 
      bkrad = 10._wp               ! light constant  

! AMMOX : oxidation of NH4 to NO2 ; light dependend
      nitrira= 0.1_wp    ! nitrification rate per day, after Yool about 0.162 per day
!             light dependency (coupled to abs_bgc, max at no light, in surface layer 0.


  !     -----------------------------------------------------------------
  !*            SET MEAN TOTAL [CA++] IN SEAWATER (MOLES/KG)
  !             (SEE BROECKER A. PENG, 1982, P. 26)
  !             ([CA++](MOLES/KG)=1.026E-2*(S/35.) AFTER
  !             CULKIN(1965), CF. BROECKER ET AL. 1982)
  !             ------------- --- -------- -- --- -----
  !
  calcon = 1.03e-2_wp

  rrrcl = salchl * 1.025_wp * bor1 * bor2
  !$ACC UPDATE DEVICE(rrrcl) ASYNC(1)

  END SUBROUTINE SET_PARAMETERS_BGC

  SUBROUTINE ini_aggregate_parameters

  ! CD parameters (formula 16)
   AJ1 = 24.00_wp
   AJ2 = 29.03_wp
   AJ3 = 14.15_wp
   BJ1 = 1.0_wp
   BJ2 = 0.871_wp
   BJ3 = 0.547_wp

   ! aggregate number distribution slope b and fractal dimension df !!!!
   ! requires: b>df+2
   agg_re_crit    = 20._wp  ! critical particle Reynolds number for limiting nr-distribution
   agg_org_dens   = 1100._wp ! detritus density - don't use orgdens to avoidnegative ws

   ! POM in HAMOCC: 122 C + 263 H + 74 O + 16 N + 1 P
   ! 122*12 + 263*1 + 74 * 16 + 16*14 + 1*31 = 3166g POM / mol organic P
   det_mol2mass   = 3166._wp ! unit: kg POM / (kmol organic P)
   dp_dust = 2.e-6_wp   ! following the classical HAMOCC parametrization
   dp_det  = 4.e-6_wp   ! not well defined
   dp_calc = 3.e-6_wp   ! following Henderiks 2008, Henderiks & Pagani 2008
   dp_opal = 20.e-6_wp  ! rough guestimate - literature search required
   stickiness_tep    = 0.19_wp
   stickiness_det    = 0.1_wp
   stickiness_opal   = 0.08_wp
   stickiness_calc   = 0.09_wp
   stickiness_dust   = 0.07_wp
   agg_df_max        = 2.4_wp
   agg_df_min        = 1.6_wp
   rho_tep           = 800._wp ! 700.-840. kg/m^3 Azetsu-Scott & Passow 2004

  END SUBROUTINE ini_aggregate_parameters

  SUBROUTINE BGC_PARAM_CONV_UNIT

   
    remido   = remido * dtb    ! 
    dyphy    = dyphy * dtb    ! 1/d -mortality rate of phytoplankton
    grazra   = grazra * dtb      ! 1/d -grazing rate 
    spemor   = spemor * dtb     ! 1/d -mortality rate of zooplankton
    gammap   = gammap * dtb     ! 1/d -exudation rate
    gammaz   = gammaz * dtb     ! 1/d -excretion rate
    sinkspeed_opal = sinkspeed_opal * dtb 
    sinkspeed_calc = sinkspeed_calc * dtb  
    sinkspeed_dust = sinkspeed_dust * dtb
    wcya = buoyancyspeed_cya *dtb  !  buoyancy speed of cya  
    sinkspeed_martin_ez = sinkspeed_martin_ez * dtb
    sinkspeed_poc = sinkspeed_poc * dtb
    sulfate_reduction = sulfate_reduction * dtb
    drempoc  = drempoc  *dtb    ! 1/d      
    dremopal = dremopal * dtb  ! 1/d      
    dremn2o  = dremn2o * dtb      ! 1/d
    dremcalc = dremcalc *dtb    ! 
    !$ACC UPDATE DEVICE(dremcalc) ASYNC(1)
    denitrification = denitrification *dtb 
    denit_sed = denit_sed *dtb    ! sediment denitrification rate
    sred_sed = sred_sed *dtb    ! sediment sulfate reduction rate
    relaxfe = relaxfe *dtb       ! relaxation time for iron to fesoly 
    disso_op = disso_op * dtb
    disso_po = disso_po * dtb
    disso_cal = disso_cal * dtb
    rh2sox  = rh2sox * dtb
    cycdec  = cycdec * dtb

    no3nh4red = no3nh4red * dtb
    no3no2red = no3no2red * dtb
    no2denit = no2denit * dtb
    anamoxra = anamoxra * dtb
    nitriox = nitriox * dtb
    nitrira = nitrira * dtb

  END SUBROUTINE

  ! ---------------------------------------------------------------------
  SUBROUTINE ini_wpoc(local_bgc_memory, ptiestw)
  ! initialize wpoc,wopal,wcal
  ! if lmartin==TRUE (mo_control_bgc, nml)
  ! wpoc increases linearly with depth below mc_depth (beleg, nml)
  ! otherwise the constant sinkspeed_poc (beleg, nml) is used
   TYPE(t_bgc_memory), POINTER :: local_bgc_memory
   REAL(wp),INTENT(in):: ptiestw(bgc_zlevs+1)

   INTEGER :: k 
   REAL(wp) :: at_mc_depth
   
   ! default case: constant sinking speed
   local_bgc_memory%wpoc(:,:)  = sinkspeed_poc
   local_bgc_memory%wopal(:,:) = sinkspeed_opal
   local_bgc_memory%wcal(:,:) = sinkspeed_calc
   local_bgc_memory%wdust(:,:) = sinkspeed_dust

   IF(i_settling==1)then
   DO k = 1,bgc_zlevs
      ! Are we at z0 (below which the Martin curve starts)?
      at_mc_depth=merge(0._wp,1._wp,ptiestw(k+1)<=mc_depth)
      ! w=w0 + a*(z-z0)
      ! z0= mc_depth
      ! a=remin_rate/b  with F(z)=F(z0)(z/zo)**(-b) 
      local_bgc_memory%wpoc(:,k) = sinkspeed_martin_ez + at_mc_depth * drempoc/mc_fac * (ptiestw(k+1) - mc_depth)
   ENDDO
   ENDIF
  END SUBROUTINE ini_wpoc
  
  SUBROUTINE ini_aquatic_tracers (local_bgc_mem, start_idx, end_idx , klevs, ibek )

    TYPE(t_bgc_memory), POINTER :: local_bgc_mem
 
    INTEGER, INTENT(in)  :: start_idx                  !< 1st REAL of model grid.
    INTEGER, INTENT(in)  :: end_idx                  !< 2nd REAL of model grid.
    INTEGER :: klevs(bgc_nproma)                  !< 3rd (vertical) REAL of model grid.
    INTEGER, INTENT(in)  :: ibek(bgc_nproma)                  !< 3rd (vertical) REAL of model grid.

    INTEGER :: j,k, m, kpke 

    REAL(wp) :: phosat, phosmed, phospac
    REAL(wp) :: oxyat, oxymed, oxypac
    REAL(wp) :: silat, silmed, silpac

    !  Initial values for aquatic (advected) ocean tracers (for new start)

    phosat  = 1.e-6_wp
    phospac = 2.5e-6_wp
    phosmed = 2.e-7_wp
    silat   = 2.e-5_wp
    silmed  = 5.e-6_wp
    silpac  = 1.e-4_wp
    oxypac  = 1.5e-4_wp
    oxyat   = 2.5e-4_wp
    oxymed  = 2.e-4_wp

 
   DO j = start_idx, end_idx
    kpke=klevs(j)
    m=ibek(j)
    DO k = 1, kpke

                local_bgc_mem%bgctra(j,k,isco212) = 2.27e-3_wp     ! [kmol/m3]
                local_bgc_mem%bgctra(j,k,ialkali) = 2.37e-3_wp
                local_bgc_mem%bgctra(j,k,iphosph) = phosmed
                local_bgc_mem%bgctra(j,k,ioxygen) = oxymed
                local_bgc_mem%bgctra(j,k,isilica) = silmed
                local_bgc_mem%bgctra(j,k,iagesc) = 0._wp
                local_bgc_mem%bgctra(j,k,igasnit)= 1e-5_wp
                local_bgc_mem%bgctra(j,k,idoc)   = 1.e-10_wp
                local_bgc_mem%bgctra(j,k,iphy)   = 1.e-8_wp
                local_bgc_mem%bgctra(j,k,izoo)   = 1.e-8_wp
                local_bgc_mem%bgctra(j,k,idet)   = 1.e-8_wp
                local_bgc_mem%bgctra(j,k,icalc)  = 0._wp
                local_bgc_mem%bgctra(j,k,ih2s)   = 0._wp
                local_bgc_mem%bgctra(j,k,idms)   = 0._wp
                local_bgc_mem%bgctra(j,k,iopal)  = 1.e-8_wp
                local_bgc_mem%bgctra(j,k,idust)  = 0._wp
                local_bgc_mem%bgctra(j,k,icya)   = 1.e-10_wp
                local_bgc_mem%bgctra(j,k,ian2o)  = 1.e-9_wp
                local_bgc_mem%hi(j,k)     = 3.e-9_wp
                !local_bgc_mem%co3(j,k)             = 0._wp       ! this good for initialisation -> 2.e-4?
                local_bgc_mem%co3(j,k)    = 2.e-4_wp       ! this good for initialisation -> 2.e-4?
                local_bgc_mem%bgctra(j,k,iano3)  = rnit*local_bgc_mem%bgctra(j,k,iphosph)
                local_bgc_mem%bgctra(j,k,iiron)  = 0.6e-9_wp

                if((m.ge.bgc_soce).and.(m.le.bgc_npac))then
                   local_bgc_mem%bgctra(j,k,iphosph) = phospac
                   local_bgc_mem%bgctra(j,k,ioxygen) = oxypac
                   local_bgc_mem%bgctra(j,k,isilica) = silpac
                   local_bgc_mem%bgctra(j,k,iano3)  = rnit*local_bgc_mem%bgctra(j,k,iphosph)
                elseif((m.ge.bgc_gin).and.(m.le.bgc_tatl))THEN
                   local_bgc_mem%bgctra(j,k,iphosph) = phosat
                   local_bgc_mem%bgctra(j,k,ioxygen) = oxyat
                   local_bgc_mem%bgctra(j,k,isilica) = silat
                   local_bgc_mem%bgctra(j,k,iano3)  = rnit*local_bgc_mem%bgctra(j,k,iphosph)
                endif

                if (l_N_cycle) then
                   local_bgc_mem%bgctra(j,k,iammo) = 1.e-2_wp*local_bgc_mem%bgctra(j,k,iano3)
                   local_bgc_mem%bgctra(j,k,iano2) = 1.e-2_wp*local_bgc_mem%bgctra(j,k,iano3)
                endif

                if(m.eq.bgc_land)then
                   local_bgc_mem%bgctra(j,k,iagesc)  = rmasko
                   local_bgc_mem%bgctra(j,k,iphosph) = rmasko
                   local_bgc_mem%bgctra(j,k,isilica) = rmasko
                   local_bgc_mem%bgctra(j,k,ioxygen) = rmasko
                   local_bgc_mem%bgctra(j,k,ialkali) = rmasko
                   local_bgc_mem%bgctra(j,k,iiron)   = rmasko
                   local_bgc_mem%bgctra(j,k,isco212) = rmasko
                   local_bgc_mem%bgctra(j,k,iano3)   = rmasko
                   local_bgc_mem%bgctra(j,k,igasnit) = rmasko
                   local_bgc_mem%bgctra(j,k,idoc)    = rmasko
                   local_bgc_mem%bgctra(j,k,iphy)    = rmasko
                   local_bgc_mem%bgctra(j,k,izoo)    = rmasko
                   local_bgc_mem%bgctra(j,k,idet)    = rmasko
                   local_bgc_mem%bgctra(j,k,icalc)   = rmasko
                   local_bgc_mem%bgctra(j,k,ih2s)    = rmasko
                   local_bgc_mem%bgctra(j,k,iopal)   = rmasko
                   local_bgc_mem%bgctra(j,k,ian2o)   = rmasko
                   local_bgc_mem%bgctra(j,k,iiron)   = rmasko
                   local_bgc_mem%bgctra(j,k,icya)    = rmasko
                   local_bgc_mem%bgctra(j,k,idms)    = rmasko
                   local_bgc_mem%hi(j,k)  = rmasko
                   local_bgc_mem%co3(j,k) = rmasko
                   if (l_N_cycle) then
                      local_bgc_mem%bgctra(j,k,iammo) = rmasko
                      local_bgc_mem%bgctra(j,k,iano2) = rmasko
                   endif
                ENDIF

       ENDDO
    ENDDO
 
 

  END SUBROUTINE ini_aquatic_tracers

  ! ---------------------------------------------------------------------

  SUBROUTINE ini_pore_water_tracers(local_bgc_mem, local_sediment_mem, start_idx,end_idx)

    TYPE(t_bgc_memory), POINTER :: local_bgc_mem
    TYPE(t_sediment_memory), POINTER :: local_sediment_mem
    
    INTEGER, INTENT(in) :: start_idx
    INTEGER, INTENT(in) :: end_idx

    INTEGER,  POINTER  :: kbo(:)   !< k-index of bottom layer (2d)    
    INTEGER :: j, k
    
    !  Initial values for sediment pore water tracers. (solid components?)
    kbo => local_bgc_mem%kbo

   DO j = start_idx, end_idx
    DO k = 1, ks
             IF(local_bgc_mem%bolay(j) > 0._wp) THEN
                local_sediment_mem%powtra(j,k,ipowaic) = local_bgc_mem%bgctra(j,kbo(j),isco212)
                local_sediment_mem%powtra(j,k,ipowaal) = local_bgc_mem%bgctra(j,kbo(j),ialkali)
                local_sediment_mem%powtra(j,k,ipowaph) = local_bgc_mem%bgctra(j,kbo(j),iphosph)
                local_sediment_mem%powtra(j,k,ipowaox) = local_bgc_mem%bgctra(j,kbo(j),ioxygen)
                local_sediment_mem%powtra(j,k,ipown2)  = 0._wp
                local_sediment_mem%powtra(j,k,ipowh2s) = 0._wp
                local_sediment_mem%powtra(j,k,ipowno3) = local_bgc_mem%bgctra(j,kbo(j),iano3)
                local_sediment_mem%powtra(j,k,ipowasi) = local_bgc_mem%bgctra(j,kbo(j),isilica)
                local_sediment_mem%powtra(j,k,ipowafe) = local_bgc_mem%bgctra(j,kbo(j),iiron)
                                
                IF (l_N_cycle) THEN
                   local_sediment_mem%powtra(j,k,ipownh4) = local_bgc_mem%bgctra(j,kbo(j),iammo)
                   local_sediment_mem%powtra(j,k,ipowno2) = local_bgc_mem%bgctra(j,kbo(j),iano2)
                ENDIF
                local_sediment_mem%sedlay(j,k,issso12) = 1.e-8_wp
                local_sediment_mem%sedlay(j,k,isssc12) = 1.e-8_wp
                local_sediment_mem%sedlay(j,k,issster) = 30._wp
                local_sediment_mem%sedlay(j,k,issssil) = 0._wp
                local_sediment_mem%sedhpl(j,k)         = local_bgc_mem%hi(j,kbo(j))
             ELSE
                local_sediment_mem%powtra(j,k,ipowno3) = rmasks   ! pore water
                local_sediment_mem%powtra(j,k,ipown2)  = rmasks
                local_sediment_mem%powtra(j,k,ipowaic) = rmasks
                local_sediment_mem%powtra(j,k,ipowaal) = rmasks
                local_sediment_mem%powtra(j,k,ipowaph) = rmasks
                local_sediment_mem%powtra(j,k,ipowh2s) = rmasks
                local_sediment_mem%powtra(j,k,ipowaox) = rmasks
                local_sediment_mem%powtra(j,k,ipowasi) = rmasks
                local_sediment_mem%powtra(j,k,ipowafe) = rmasks
                IF (l_N_cycle) THEN
                   local_sediment_mem%powtra(j,k,ipownh4) = rmasks
                   local_sediment_mem%powtra(j,k,ipowno2) = rmasks
                ENDIF
                local_sediment_mem%sedlay(j,k,issso12) = rmasks   ! solid sediment
                local_sediment_mem%sedlay(j,k,isssc12) = rmasks
                local_sediment_mem%sedlay(j,k,issssil) = rmasks
                local_sediment_mem%sedlay(j,k,issster) = rmasks
                local_sediment_mem%sedhpl(j,k)         = rmasks
             ENDIF

       ENDDO
    ENDDO
 
 

  END SUBROUTINE ini_pore_water_tracers

  ! ---------------------------------------------------------------------
  SUBROUTINE ini_continental_carbon_input(totarea)

    REAL(wp), INTENT(in):: totarea

    calcinp = deltacalc * dtbgc / totarea
    orginp  = deltaorg  * dtbgc / totarea
    silinp  = deltasil  * dtbgc / totarea

  END SUBROUTINE ini_continental_carbon_input

  SUBROUTINE ini_atmospheric_concentrations(local_bgc_mem)

    TYPE(t_bgc_memory), POINTER :: local_bgc_mem
    !  atmospheric concentrations (overwritten from restart file)

    local_bgc_mem%atm(:,iatmco2) = atm_co2
    local_bgc_mem%atm(:,iatmo2)  = atm_o2
    local_bgc_mem%atm(:,iatmn2)  = atm_n2

  END SUBROUTINE ini_atmospheric_concentrations

END MODULE 
