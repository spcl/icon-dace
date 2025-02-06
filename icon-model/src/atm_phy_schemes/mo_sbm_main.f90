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

! Description:
!
! this is the spectral-bin microphysics scheme based on the hebrew university
! cloud model (hucm), originally formulated and coded by alexander khain
!
! the wrf bin microphysics scheme (fast sbm or fsbm) solves equations for four
! size distribution FUNCTIONs: aerosols, drop (including rain drops), snow and
! graupel/hail (from which mass mixing ratio qna, qc, qr, qs, qg/qh and
! their number concentrations are calculated).

! the scheme is generally written in cgs units. in the updated scheme (fsbm-2)
! the users can choose either graupel or hail to describe dense particles
! (see the 'hail_opt' switch). by default, the 'hail_opt = 1' is used.
! hail particles have larger terminal velocity than graupel per mass bin.
! 'hail_opt' is recommended to be used in simulations of continental clouds
! systems. the graupel option may lead to better results in simulations of
! maritime convection.

! the aerosol spectrum in fsbm-2 is approximated by 3-lognormal size distribution
! representing smallest aerosols (nucleation mode), intermediate-size
! (accumulation mode) and largest aerosols (coarse mode). the bc/ic for aerosols
! ,as well as aerosols vertical distribution profile -- are set from within the
! fsbm-2 scheme (see the 'domain_id' parameter). the domain_id forces bc to be applied
! for the parent domain only.

! the user can set the liquid water content threshold (lwc) in which rimed snow
! is being transferred to hail/graupel (see 'alcr' parameter).
! the default value is alcr = 0.5 [g/m3]. increasing this value will result
! in an increase of snow mass content, and a decrease in hail/graupel mass
! contents.

! we thank and acknowledge contribution from jiwen fan (pnnl), alexander rhyzkov
! (cimms/nssl), jeffery snyder (cimms/nssl), jimy dudhia (ncar) and dave gill! (ncar).

! useful references:
! -------------------
!   khain, a. p., and i. sednev, 1996: simulation of precipitation formation in
! the eastern mediterranean coastal zone using a spectral microphysics cloud
! ensemble model. atmospheric research, 43: 77-110;
!   khain, a. p., a. pokrovsky and m. pinsky, a. seifert, and v. phillips, 2004:
! effects of atmospheric aerosols on deep convective clouds as seen from
! simulations using a spectral microphysics mixed-phase cumulus cloud model
! part 1: model description. j. atmos. sci 61, 2963-2982);
!   khain a. p. and m. pinsky, 2018: physical processes in clouds and cloud
! modeling. cambridge university press. 642 pp
!   shpund, j., a. khain, and d. rosenfeld, 2019: effects of sea spray on the
! dynamics and microphysics of an idealized tropical cyclone. j. atmos. sci., 0,
! https://doi.org/10.1175/jas-d-18-0270.1 (a preliminary description of the
! updated fsbm-2 scheme)

! when using the fsbm-2 version please cite:
! -------------------------------------------
! shpund, j., khain, a., lynn, b., fan, j., han, b., ryzhkov, a., snyder, j.,
! dudhia, j. and gill, d., 2019. simulating a mesoscale convective system using wrf
! with a new spectral bin microphysics: 1: hail vs graupel.
! journal of geophysical research: atmospheres.
! note by jiwen fan
! (1) the main SUBROUTINE is fast_sbm where all the microphysics processes are
!     called
! (2) for aerosol setup, seach "aerosol setup", where one can set up aerosol sd,
! composition information (molecular weight, ions, and density). for sd, there
! is a choice for a lognormal distribution, or read from an observed sd.
! (3) my postdoc yuwei zhang has added cloud related diagnostics (mainly process
! rates) and added an option to read in the observed sd. observed sd data should be
! processed following! a format of the file "ccn_size_33bin.dat" which is in
! size (cm), dn (# cm-3), and dndlogd for 33bins

MODULE mo_sbm_main

USE mo_kind,                 ONLY: wp, dp
USE mo_exception,            ONLY: finish, message, txt => message_text
!USE mo_output_event_types,   ONLY: t_event_step
!USE mo_run_config,           ONLY: iqbin, iqb_i, iqb_e, iqb_s

USE mo_sbm_util,             ONLY:                                                  &
    ro1bl, xl, vr1, use_cloud_base_nuc, t_nucl_drop_min, isign_3point, isign_ko_1,  &
    isign_ko_2, coeff_remaping, ratio_icew_min, rw_pw_min, rw_pw_ri_pi_min,         &
    ventpl_max, g_lim, kp_flux_max, accr_cld_msink_b, accr_cld_nsink_b,             &
    auto_cld_msink_b, auto_cld_nsink_b, chucm, col, cwll, nkr, icemax, jmax,        &
    ibreakup, ima, jbreak, krmin_breakup, pkij, qkj, selfc_rain_nchng_b,            &
    lh_ce_1, ccn_reg, ce_af, ce_bf, cldnucl_af, cldnucl_bf,                         &
    conserv, del_ccnreg, del_ce_sum, del_cldnucl_sum, del_ds_sum,                   &
    dt_coll, ions, krdrop, mwaero, ncoll, ncond, nkr_aerosol, p_ff1i01, p_ff1i33,   &
    p_ff8i01, p_ff8i33, rccn, rlec, ro_solute, xl_mg,                               &
    kernals_ks


USE mo_physical_constants,   ONLY: cpd, cvd
                                   ! [J/K/kg] specific heat at constant pressure
                                   ! [J/K/kg] specific heat at constant volume
!USE mo_satad,                ONLY: qsat_rho

  PRIVATE
  PUBLIC warm_sbm
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_sbm_main'
  REAL(KIND=wp) :: r_v = 461.51_wp

  CONTAINS

  SUBROUTINE warm_sbm (dt                &!in:    dt
                      ,dz8w              &!in:    vertical layer thickness
                      ,rho_phy           &!in:    density
                      ,p_phy             &!in:    pressure
                      ,pi_phy            &!in:    exner
                      ,w                 &!in:    updraft velocity
                      ,qv_old            &!in     cloud water before advection
                      ,th_phy            &!inout: theta. Check how to update prognostic theta_v ???
                      ,qv                &!inout: specific humidity
                      ,chem_new          &!inout: 99 mass bins
                      ,rainncv           &!inout: 1 time step precip. (mm/sec): input: 0, output can go further to the model
                      ,qc                &!inout: cloud water:                   input: 0
                      ,qr                &!inout: rain water:                    input: 0
                      ,qnc               &!inout: cloud water concentration:     input: 0
                      ,qnr               &!inout: rain water concentration:      input: 0
                      ,qna               &!inout: ccn concentration:             input: 0
                      ,qna_nucl          &!inout: nucleated ccn concentration:   input: 0
                      ,lh_rate           &!inout: rate 1:                        input: 0, output can    go further to the model
                      ,ce_rate           &!inout: rate 2:                        input: 0, output can    go further to the model
                      ,cldnucl_rate      &!inout: rate 3:                        input: 0, output can    go further to the model
                      ,its,ite, kts,kte  &!in:    subdomain indeces
                      ,diag_satur_ba,diag_satur_aa,diag_satur_am,temp_old &
                      ,temp_new,reff,reffc,reffr, diag_supsat_out)

  IMPLICIT NONE

  INTEGER :: kr,ikl
  INTEGER,INTENT(IN) :: its,ite,kts,kte
  REAL(KIND=wp), INTENT(IN)     :: dt
  REAL(KIND=wp), DIMENSION(:,:), INTENT(IN) ::   w
  REAL(KIND=wp), DIMENSION(:,:,:),INTENT(INOUT)   :: chem_new
  REAL(KIND=wp), DIMENSION(its:ite, kts:kte), INTENT(INOUT) :: & 
   qc,       &
   qnc,      &
   qr,       &
   qnr,      &
   qna,      &
   qna_nucl, &
   reff,reffc,reffr

  REAL(KIND=wp), DIMENSION(:,:), INTENT(INOUT) :: lh_rate,ce_rate,cldnucl_rate
  REAL(KIND=wp), DIMENSION(:,:), INTENT(IN) ::  qv_old, temp_old, temp_new
  REAL(KIND=wp), DIMENSION(its:ite, kts:kte) ::  qv_oldtmp
  REAL(KIND=wp), INTENT(INOUT), DIMENSION(:) :: rainncv
  REAL(KIND=wp), INTENT(IN), DIMENSION(:,:) :: dz8w,p_phy,pi_phy,rho_phy
  REAL(KIND=wp), INTENT(INOUT),  DIMENSION(:,:) :: th_phy, qv,diag_satur_ba,diag_satur_aa,diag_satur_am,diag_supsat_out
  REAL(KIND=wp),  DIMENSION(its:ite, kts:kte) :: t_new,t_old,zcgs,rhocgs,pcgs
  REAL(KIND=wp),  DIMENSION(1:33) :: fccn, fccn_nucl
  REAL(KIND=wp) :: mom_2c,mom_3c,mom_2r,mom_3r,mom_2,mom_3,supsat_out
  INTEGER :: i,k !,n_chem
  REAL(KIND=wp) ::  &
!     &        sup2_old, &
      &        aa1_my,bb1_my,aa2_my,bb2_my, &
      &        dtcond,dtnew,dtcoll, &
      &        a1_myn, bb1_myn, a2_myn, bb2_myn
  DATA a1_myn, bb1_myn, a2_myn, bb2_myn  &
      &      /2.53,5.42,3.41e1,6.13/
  DATA aa1_my,bb1_my,aa2_my,bb2_my/2.53e12,5.42e3,3.41e13,6.13e3/
  REAL(KIND=wp),DIMENSION (nkr) :: ff1in,ff1r


  REAL(KIND=wp) :: del1nr,del2nr,del12r,del12rd,es1n,es2n,ew1n,ew1pn
  REAL(KIND=wp) :: delsup1,delsup2,deldiv1,deldiv2
  REAL(KIND=wp) :: tt,qq,tta,qqa,pp
  REAL(KIND=wp) :: div1,div2,div3,div4,del1in,del2in,del1ad,del2ad
  REAL(KIND=wp) :: del_bb,del_bbn,del_bbr
  INTEGER krr,i_start,i_end
  REAL(KIND=wp) :: fmax1
  INTEGER isym1,isym2(icemax),isym3,isym4,isym5
  INTEGER diffu
  REAL(KIND=wp) :: deltaw
  REAL(KIND=wp) :: zcgs_z(kts:kte),pcgs_z(kts:kte),rhocgs_z(kts:kte),ffx_z(kts:kte,nkr)
  REAL(KIND=wp) :: z_full
  REAL(KIND=wp) :: vrx(kts:kte,nkr)
  REAL(KIND=wp) :: vr1_z(nkr,kts:kte), factor_p, vr1_z3d(nkr,its:ite,kts:kte)
  REAL(KIND=wp), PARAMETER :: pi=3.14159265359
  REAL(KIND=wp) ::    w_stag_my
  ! ... polar-hucm
  INTEGER :: is_this_cloudbase
  INTEGER :: sat_method
  REAL(KIND=wp) :: satur_diag,ew_tmp 
  LOGICAL :: latheatfac1, latheatfac2, latheatfac3, latheatfac4

  latheatfac1=.FALSE.
  latheatfac2=.TRUE.   
  latheatfac3=.TRUE.   !paper
  latheatfac4=.FALSE.

  sat_method=3
 !n_chem=iqb_e

  ncond = 3
  ncoll = 1
  dtcond = dt/REAL(ncond) ! option for long relaxation: dtcond = 200.0*dt/real(ncond)
  dtcoll = dt/REAL(ncoll)
  dt_coll = dtcoll

  del_bb=bb2_my-bb1_my
  del_bbn=bb2_myn-bb1_myn
  del_bbr=bb1_myn/del_bbn
  i_start=its
  i_end=ite

  DO i = i_start,i_end
    DO k = kts,kte
      qv_oldtmp(i,k)=qv_old(i,k)
    END DO
  END DO
  DO k = kts,kte
    DO i = i_start,i_end
      t_new(i,k)=temp_new(i,k)
      t_old(i,k) = temp_old(i,k)
    END DO
  END DO

  DO i = its,ite
    DO k = kts,kte

      rhocgs(i,k)=rho_phy(i,k)*0.001
    ! ... drops
      krr=0
      DO kr=p_ff1i01,p_ff1i33
        krr=krr+1
        chem_new(i,k,kr)=chem_new(i,k,kr)*rhocgs(i,k)/col/xl(krr)/xl(krr)/3.0 ! input from model: kg/kg. inside: #/(g*cm^3)
      END DO
    ! ... aerosols
      krr=0
      DO kr=p_ff8i01,p_ff8i33
        krr=krr+1
        chem_new(i,k,kr) = chem_new(i,k,kr)*rhocgs(i,k)/1000.0   !input from model: #/kg. inside: #/cm^3
      END DO
    ! ... nucleated aerosols - not implemented here
    END DO
  END DO

  DO i = i_start,i_end
    z_full=0.
    DO k = kte,kts,-1 ! was do k = kts,kte before reversing z axis
      pcgs(i,k)=p_phy(i,k)*10.0_wp
      rhocgs(i,k)=rho_phy(i,k)*0.001_wp
      zcgs(i,k)=z_full+0.5*dz8w(i,k)*100.0_wp
      z_full=z_full+dz8w(i,k)*100.0_wp
    END DO
  END DO

  DO k = kts,kte
    DO i = its,ite
      ! ... correcting look-up-table terminal velocities (stronger with height):
      factor_p = dsqrt(1.0d6/pcgs(i,k))
      vr1_z(1:nkr,k) =  vr1(1:nkr)*factor_p
      vr1_z3d(1:nkr,i,k) = vr1(1:nkr)*factor_p

      ! ... droplet / drops
      krr = 0
      DO kr = p_ff1i01,p_ff1i33
        krr = krr + 1
        ff1r(krr) = chem_new(i,k,kr)
        IF (ff1r(krr) < 0.0)ff1r(krr) = 0.0
      END DO
      ! ... ccn
      krr = 0
      DO kr=p_ff8i01,p_ff8i33
        krr = krr + 1
        fccn(krr) = chem_new(i,k,kr)
        IF (fccn(krr) < 0.0)fccn(krr) = 0.0
      END DO           
      ! ... nucleated ccn - not implemented here
      lh_ce_1 = 0.0;
      ce_bf = 0.0; ce_af = 0.0; cldnucl_af = 0.0; cldnucl_bf = 0.0; 
      del_cldnucl_sum = 0.0; del_ce_sum = 0.0; del_ds_sum=0.0;
            
      auto_cld_nsink_b = 0.0; auto_cld_msink_b = 0.0;  
      accr_cld_nsink_b = 0.0; accr_cld_msink_b = 0.0;
      selfc_rain_nchng_b = 0.0
      
     ! +---------------------------------------------+
     ! neucliation, condensation, collisions
     ! +---------------------------------------------+
     satur_diag=0.0
     CALL satcalc_diag(t_old(i,k),qv_oldtmp(i,k),pcgs(i,k),aa1_my,bb1_my,satur_diag,sat_method,rho_phy(i,k))
     diag_satur_ba(i,k)=satur_diag
     CALL satcalc_diag(t_new(i,k),qv(i,k),pcgs(i,k),aa1_my,bb1_my,satur_diag,sat_method,rho_phy(i,k))
     diag_satur_aa(i,k)=satur_diag

     ccn_reg = 0.0
     del_ccnreg = 0.0  
     tt=t_old(i,k)  !temperature before advection
     qq=qv_oldtmp(i,k) !specific humidity before advection kg/kg
     IF (qq.LE.0.0) qq = 1.d-10
     pp=pcgs(i,k)
     tta=t_new(i,k) !temperature after advection
     qqa=qv(i,k)    !specific humidity after advection

     IF (qqa.LE.0) CALL message(modname,"warning: fast sbm, qqa < 0")
       IF (qqa.LE.0) PRINT*,'i,k,told,tnew,qqa = ',i,k,tt,tta,qqa
       IF (qqa.LE.0) qqa = 1.0d-10

       IF (sat_method == 1) THEN
         es1n = aa1_my*EXP(-bb1_my/tt)
         es2n = aa2_my*EXP(-bb2_my/tt)
       ELSE
         es1n=10.0*610.78*EXP(17.269*(tt-273.15)/(tt-35.86))
         es2n=10.0*610.78*EXP(21.875*(tt-273.15)/(tt-7.66))
       END IF
       IF ((sat_method == 1) .OR. (sat_method == 2)) THEN
         ew1n=qq*pp/(0.622+0.378*qq)
       ELSE
         ew1n=10.0*qq*rho_phy(i,k)*r_v*tt !icon: qq[kg/kg], rho[kg/m3]
       END IF

       div1=ew1n/es1n
       del1in=ew1n/es1n-1. !supersaturation over water before advection
       div2=ew1n/es2n      !supersaturation over ice before advection
       del2in=ew1n/es2n-1.

       IF (sat_method == 1) THEN
         es1n=aa1_my*EXP(-bb1_my/tta)
         es2n=aa2_my*EXP(-bb2_my/tta)
       ELSE
         es1n=10.0*610.78*EXP(17.269*(tta-273.15)/(tta-35.86))
         es2n=10.0*610.78*EXP(21.875*(tta-273.15)/(tta-7.66))
       END IF
       IF ((sat_method == 1) .OR. (sat_method == 2)) THEN
         ew1n=qqa*pp/(0.622+0.378*qqa)
       ELSE
         ew1n=10.0*qqa*rho_phy(i,k)*r_v*tta !icon
       END IF

       div3=ew1n/es1n
       del1ad=ew1n/es1n-1. !supersaturation over water after advection

       div4=ew1n/es2n
       del2ad=ew1n/es2n-1. !supersaturation over ice after advection
!      sup2_old=del2in

       IF (del1ad > 0.0 .OR. (sum(ff1r)) > 0.0) THEN
         ! previously there was an option: "call relaxation_time" here
         delsup1=(del1ad-del1in)/ncond
         delsup2=(del2ad-del2in)/ncond
         deldiv1=(div3-div1)/ncond
         deldiv2=(div4-div2)/ncond
                
         diffu=1
         IF (div1.EQ.div3) diffu = 0
         IF (div2.EQ.div4) diffu = 0 

         dtnew = 0.0
         DO ikl=1,ncond
           dtcond = min(dt-dtnew,dtcond)
           !option for long relaxation: dtcond = min(200.0*dt-dtnew,dtcond)
           dtnew = dtnew + dtcond

           IF (diffu == 1)THEN
             IF (diffu == 1)THEN
               del1in = del1in+delsup1
               del2in = del2in+delsup2
               div1 = div1+deldiv1
               div2 = div2+deldiv2
             END IF
             IF ((div1.GE.div2.and.tt.LE.272.0) .OR. (tt.LT.223.0))THEN
               diffu = 0
             END IF

             IF (diffu == 1)THEN
               del1nr=a1_myn*(100.*div1)           !move to method=1
               del2nr=a2_myn*(100.*div2)           !move to method=1
               IF (del2nr.EQ.0)PRINT*,'ikl = ',ikl
               IF (del2nr.EQ.0)PRINT*,'div1,div2 = ',div1,div2
               IF (del2nr.EQ.0)PRINT*,'i,k = ',i,k
               IF (del2nr.EQ.0)CALL finish(TRIM(modname),"fatal error in module_mp_fast_sbm (del2nr.EQ.0) , model stop ")
               del12r=del1nr/del2nr                    !move to method=1
               del12rd=del12r**del_bbr                 !move to method=1
               ew1pn=aa1_my*100.*div1*del12rd/100.     !move to method=1
               IF (sat_method == 1) THEN
                 tt=-del_bb/DLOG(del12r) !tt for given supersaturation
                 qq=0.622*ew1pn/(pp-0.378*ew1pn) !qq for given supersaturation
               ELSE
                 CALL temp_from_sat(tt,del1in,del2in) !,qq,rho_phy(i,k)) !pp) !del1in,del2in: supsat over water and ice
               END IF
               IF (sat_method == 2) THEN
                 ew_tmp=(del1in+1.0)*10.0*610.78*EXP(17.269*(tt-273.15)/(tt-35.86))
                 qq=0.622*ew_tmp/(pp-0.378*ew_tmp)
               ELSE IF (sat_method == 3) THEN
                 ew_tmp=(del1in+1.0)*10.0*610.78*EXP(17.269*(tt-273.15)/(tt-35.86))
                 qq=ew_tmp/(rho_phy(i,k)*r_v*tt*10.0) !icon
               END IF

               IF (del1in .GT. 0.0 .OR. del2in .GT. 0.0)THEN
                 ! +------------------------------------------+
                 ! droplet nucleation :
                 ! +------------------------------------------+
                 ff1in(:) = ff1r(:)
                 cldnucl_bf = 3.0*col*( sum(ff1in*(xl**2.0)) )/rhocgs(i,k)
                 is_this_cloudbase = 0
                 w_stag_my = 0.0d0                 

                 CALL jernucl01_ks(ff1in,fccn,fccn_nucl         &
                                ,xl,tt                          &
                                ,rhocgs(i,k),pcgs(i,k)          &
                                ,del1in                         & !,del2in &
                                ,col                            &
                                ,rccn,nkr,nkr_aerosol           &
                                ,w_stag_my,is_this_cloudbase,ro_solute,ions,mwaero)

                 cldnucl_af = 3.0*col*( sum(ff1in*(xl**2.0)) )/rhocgs(i,k)
                 del_cldnucl_sum =  del_cldnucl_sum + (cldnucl_af - cldnucl_bf)
                 ff1r(:) = ff1in(:)
               END IF

               fmax1=0.0
               DO kr=1,nkr
                 ff1in(kr)=ff1r(kr)
                 fmax1=MAX(ff1r(kr),fmax1)
               END DO
               isym1 = 0
               IF (fmax1 > 0)isym1 = 1

               ce_bf = 3.0*col*( sum(ff1r*(xl**2.0)) )/rhocgs(i,k)                          
               IF ((tt>273.15) .OR. (tt<273.15-0.187)) THEN
               ! ... only warm phase - diffusional growth
                 supsat_out=0.0
                 CALL onecond1(tt,qq,pp,rhocgs(i,k) &
                              ,vr1_z(:,k)           &
                              ,del1in,del2in,div1,div2 &
                              ,ff1r,ff1in,xl,rlec,ro1bl &
                              ,aa1_my,bb1_my,aa2_my,bb2_my &
                              ,col,dtcond,icemax,nkr,isym1 &
                              ,isym2,isym3,isym4,isym5,i,k,w(i,k),ccn_reg &
                              ,latheatfac1,latheatfac2,latheatfac3,latheatfac4,supsat_out,sat_method,rho_phy(i,k))
               END IF
               ce_af = 3.0*col*( sum(ff1r*(xl**2.0)) )/rhocgs(i,k)               
               del_ce_sum = del_ce_sum + (ce_af - ce_bf)
             END IF ! diff.NE.0
           END IF ! diffu.NE.0
         END DO ! ncond - end of ncond loop

         ! switched off ccn_regeneration
        !n_reg_ccn_bf = col*sum(fccn)
        !n_reg_ccn_af = col*sum(fccn)

         ! +----------------------------------+
         ! collision-coallescnce
         ! +----------------------------------+

         DO ikl = 1,ncoll
           CALL coal_bott_new_warm (ff1r,tt,pp, dt_coll, krdrop)
         END DO ! ncoll - end of ncoll loop
                                   
         t_new(i,k) = tt
         qv(i,k) = qq !kg/kg

       END IF 
    
       !calculate supersaturation after microphysics:
       CALL satcalc_diag(t_new(i,k),qv(i,k),pcgs(i,k),aa1_my,bb1_my,satur_diag,sat_method,rho_phy(i,k))
       diag_satur_am(i,k)=satur_diag
       diag_supsat_out(i,k)=supsat_out

       ! ... process rate (integrated)
       lh_rate(i,k) = lh_rate(i,k) +  lh_ce_1/dt
       ce_rate(i,k) = ce_rate(i,k) +    del_ce_sum/dt
       cldnucl_rate(i,k) = cldnucl_rate(i,k) + del_cldnucl_sum/dt

       ! update temperature at the end of mp
       th_phy(i,k) = t_new(i,k)/pi_phy(i,k)

       ! ... drops
       krr = 0
       DO kr = p_ff1i01,p_ff1i33
         krr = krr+1
         chem_new(i,k,kr) = ff1r(krr)
       END DO
       ! ... ccn
       krr = 0
       DO kr=p_ff8i01,p_ff8i33
         krr=krr+1
         chem_new(i,k,kr)=fccn(krr)
       END DO
       ! ... nucleated ccn - not implemented here
     END DO
   END DO

   ! +-----------------------------+
   ! hydrometeor sedimentation
   ! +-----------------------------+ 
   DO i = its,ite
   ! ... drops ...
     DO k = kts,kte
       rhocgs_z(k)=rhocgs(i,k)
       pcgs_z(k)=pcgs(i,k)
       zcgs_z(k)=zcgs(i,k)
       vrx(k,:)=vr1_z3d(:,i,k)
       krr = 0
       DO kr=p_ff1i01,p_ff1i33
         krr=krr+1
         ffx_z(k,krr)=chem_new(i,k,kr)/rhocgs(i,k)
       END DO
     END DO
     CALL falfluxhucm_z(ffx_z,vrx,rhocgs_z,zcgs_z,dt,kts,kte,nkr)
     DO k = kts,kte
       krr = 0
       DO kr=p_ff1i01,p_ff1i33
         krr=krr+1
         chem_new(i,k,kr)=ffx_z(k,krr)*rhocgs(i,k)
       END DO
     END DO
   END DO

   ! ... output block
   DO k = kts,kte
     DO i = its,ite
       qc(i,k) = 0.0
       qr(i,k) = 0.0
       qnc(i,k) = 0.0
       qnr(i,k) = 0.0
       qna(i,k) = 0.0
       qna_nucl(i,k) = 0.0
       reff(i,k) = 0.0
       reffc(i,k) = 0.0
       reffr(i,k) = 0.0
       mom_2c = 0.0
       mom_3c = 0.0
       mom_2r = 0.0
       mom_3r = 0.0
       mom_2 = 0.0
       mom_3 = 0.0

       tt= th_phy(i,k)*pi_phy(i,k)

       ! ... drop output
       krr = 0
       DO kr = p_ff1i01,p_ff1i33
         krr=krr+1
         IF (krr < krdrop)THEN
           qc(i,k) = qc(i,k) &
           + (1./rhocgs(i,k))*col*chem_new(i,k,kr)*xl(krr)*xl(krr)*3    ! [qc]=kg/kg, [chem_new1-33]=#/(gr*cm^3)
           qnc(i,k) = qnc(i,k) &
           + col*chem_new(i,k,kr)*xl(krr)*3.0/rhocgs(i,k)*1000.0        ! [qnc]=#/kg, [chem_new1-33]=#/(gr*cm^3)

           mom_2c=mom_2c+3.0*xl(krr)*chem_new(i,k,kr)*(3.0*xl(krr)/(4.0*pi))**(2.0/3.0) 
           mom_3c=mom_3c+3.0*xl(krr)*chem_new(i,k,kr)*(3.0*xl(krr)/(4.0*pi)) 
         ELSE
           qr(i,k) = qr(i,k) &
           + (1./rhocgs(i,k))*col*chem_new(i,k,kr)*xl(krr)*xl(krr)*3.0  ! [qr]=kg/kg, [chem_new1-33]=#/(gr*cm^3)
           qnr(i,k) = qnr(i,k) &
           + col*chem_new(i,k,kr)*xl(krr)*3.0/rhocgs(i,k)*1000.0        ! [qnr]=#/kg, [chem_new1-33]=#/(gr*cm^3)

           mom_2r=mom_2r+3.0*xl(krr)*chem_new(i,k,kr)*(3.0*xl(krr)/(4.0*pi))**(2.0/3.0) 
           mom_3r=mom_3r+3.0*xl(krr)*chem_new(i,k,kr)*(3.0*xl(krr)/(4.0*pi)) 
         END IF
         mom_2=mom_2+3.0*xl(krr)*chem_new(i,k,kr)*(3.0*xl(krr)/(4.0*pi))**(2.0/3.0) 
         mom_3=mom_3+3.0*xl(krr)*chem_new(i,k,kr)*(3.0*xl(krr)/(4.0*pi)) 
       END DO

       IF (qc(i,k) > 1.0e-6 .and. mom_2c > 0.0) THEN
         reffc(i,k)=(mom_3c/mom_2c)*1.0e4
       END IF
       IF (qr(i,k) > 1.0e-6 .and. mom_2r > 0.0) THEN
         reffr(i,k)=(mom_3r/mom_2r)*1.0e4
       END IF
       IF (qc(i,k) + qr(i,k) > 1.0e-6 .and. mom_2 > 0.0) THEN
         reff(i,k)=(mom_3/mom_2)*1.0e4
       END IF

       ! ... aerosols output 
       krr = 0
       DO kr = p_ff8i01,p_ff8i33
         krr = krr + 1
         qna(i,k) = qna(i,k) + col*chem_new(i,k,kr)/rhocgs(i,k)*1000.0              ! [qna]=#/kg, [chem_new34-66]=#/(cm^3)
       END DO

       ! ... nucleated aerosols output - not implemented here
       qna_nucl(i,k) = 0.0
     END DO
   END DO

   DO i = its,ite
     rainncv(i) = 0.0
     krr = 0
! avoid division by zero
!$NEC novector
     DO kr=p_ff1i01,p_ff1i33
       krr=krr+1
       deltaw = vr1_z3d(krr,i,kte-1)
       rainncv(i) = rainncv(i) &
       +10.0*(3./ro1bl(krr))*col*deltaw* & !+10.0*(3./ro1bl(krr))*col*dt*deltaw* &
       chem_new(i,kte-1,kr)*xl(krr)*xl(krr)
     END DO
   END DO

   IF (conserv)THEN
     DO i = its,ite
       DO k = kts,kte        
         rhocgs(i,k)=rho_phy(i,k)*0.001
         ! ... drops  
         krr=0
         DO kr=p_ff1i01,p_ff1i33
           krr=krr+1
           chem_new(i,k,kr)=chem_new(i,k,kr)/rhocgs(i,k)*col*xl(krr)*xl(krr)*3.0
           IF (qc(i,k)+qr(i,k).LT.1.e-13)chem_new(i,k,kr)=0.0
         END DO
         ! ... ccn
         krr=0
         DO kr=p_ff8i01,p_ff8i33
           krr=krr+1
           chem_new(i,k,kr)=chem_new(i,k,kr)/rhocgs(i,k)*1000.0
         END DO
         ! ... nucleated ccn - not implemented here
         END DO
       END DO
     END IF

     RETURN
   END SUBROUTINE warm_sbm

   SUBROUTINE satcalc_diag(t,q,p,a,b,satur_diag,sat_method,rho) ! saturation calculation
     IMPLICIT NONE
     REAL(KIND=wp), INTENT(IN)   :: t,p,a,b
     REAL(KIND=wp),INTENT(IN) :: rho
     INTEGER, INTENT(IN)   :: sat_method 
     REAL(KIND=wp), INTENT(INOUT):: q
     REAL(KIND=wp), INTENT(INOUT):: satur_diag
     REAL(KIND=wp) :: es1n, ew1n

     IF (q.LE.0.0) q = 1.d-10
     IF (sat_method == 1) THEN
       es1n = a*EXP(-b/t)
     ELSE
       es1n=10.0*610.78*EXP(17.269*(t-273.15)/(t-35.86))
     END IF
     IF ((sat_method == 1) .OR. (sat_method == 2)) THEN
       ew1n=q*p/(0.622+0.378*q)
     ELSE
       ew1n=10.0*q*rho*r_v*t !icon
     END IF

     satur_diag=ew1n/es1n-1.
   END SUBROUTINE satcalc_diag

   SUBROUTINE temp_from_sat(temp,sat_w,sat_i)
     IMPLICIT NONE
     REAL(KIND=wp), INTENT(INOUT) :: sat_w,sat_i
     REAL(KIND=wp), INTENT(INOUT):: temp
     REAL(KIND=wp) :: b2w,b2i,b3,b4w,b4i, & !b1
       & a1_sat,a2_sat,a3_sat,a4_sat,a5_sat,a6_sat
   
!    b1=10.0*610.78
     b2w=17.269
     b2i=21.875
     b3=273.15
     b4w=35.86
     b4i=7.66

     a1_sat=b2w/b2i
     a2_sat=(1.0_wp/b2i)*LOG((sat_w+1.0_wp)/(sat_i+1.0_wp))
     a3_sat=a1_sat*b4i-b4w
     a4_sat=b3*(1.0-a1_sat)-a2_sat*(b4i+b4w)-a3_sat
     a5_sat=1.0-a1_sat-a2_sat
     a6_sat=b3*a3_sat+b4i*b4w*a2_sat
   
     IF (a5_sat .EQ. 0.0) THEN
       temp=-a6_sat/a4_sat
     ELSE
       IF (a4_sat**2.0+4.0*a5_sat*a6_sat .LT. 0.0) THEN
         PRINT*,'error: ',sat_w,sat_i
       END IF
       temp=(a4_sat+(a4_sat**2.0+4.0*a5_sat*a6_sat)**0.5)/(2.0*a5_sat)
     END IF
   END SUBROUTINE temp_from_sat

   SUBROUTINE falfluxhucm_z(chem_new,vr1,rhocgs,zcgs,dt,kts,kte,nkr)
     IMPLICIT NONE
     INTEGER,INTENT(IN) :: kts,kte,nkr
     REAL(KIND=wp),INTENT(INOUT) :: chem_new(:,:)
     REAL(KIND=wp),INTENT(IN) :: rhocgs(:),zcgs(:),vr1(:,:),dt
     INTEGER :: k,kr
     REAL(KIND=wp) :: tfall,dtfall,vfall(kte),dwflux(kte)
     INTEGER :: ifall,n,nsub
     ! falling fluxes for each kind of cloud particles: c.g.s. unit
     ! adapted from gsfc code for hucm
     ! the flux at k=1 is assumed to be the ground so flux(1) is the
     ! flux into the ground. dwflux(1) is at the lowest half level where
     ! q(1) etc are defined. the formula for flux(1) uses q(1) etc which
     ! is actually half a grid level above it. this is what is meant by
     ! an upstream method. upstream in this case is above because the
     ! velocity is downwards.
     ! use upstream method (vfall is positive)

     DO kr=1,nkr
       ifall=0
       DO k = kts,kte
         IF (chem_new(k,kr).GE.1.e-20) ifall=1   ! if there is some mass to fall
       END DO
       IF (ifall.EQ.1)THEN
         tfall=1.e10
         DO k=kts,kte
           ! [ks] vfall(k) = vr1(k,kr)*sqrt(1.e6/pcgs(k))
           vfall(k) = vr1(k,kr) ! ... [ks] : the pressure effect is taken into account at the beggining of the calculations
           tfall=MIN(tfall,zcgs(k)/(vfall(k)+1.e-20))   ! set tfall=z/v
         END DO
         IF (tfall.GE.1.e10) stop
         nsub=(int(2.0*dt/tfall)+1)
         dtfall=dt/nsub         ! dt used in sedimintation, which is smaller than model dt to obey CFL criterion

         DO n=1,nsub                                      ! loop over sedimintation sub steps
           DO k=kte,kts+1,-1 !loop from 2 down to 65, each level -(lower-upper)/upper. Used to be do k=kts,kte-1
            dwflux(k)=-(rhocgs(k)*vfall(k)*chem_new(k,kr)-&               ! flux difference (divergence) at level k
                        rhocgs(k-1)*vfall(k-1)*chem_new(k-1,kr))/&
                         (rhocgs(k)*(zcgs(k-1)-zcgs(k)))
           END DO
           ! no z above top, so use the same deltaz
           dwflux(kts)=-(rhocgs(kts)*vfall(kts)*chem_new(kts,kr))/&         ! flux difference (divergence) at top
                       (rhocgs(kts)*(zcgs(kts)-zcgs(kts+1)))
           DO k=kts,kte
             chem_new(k,kr)=chem_new(k,kr)+dwflux(k)*dtfall     !if the droplets fall into some level, add them to the psd there
           END DO
         END DO
       END IF
     END DO

     RETURN
   END SUBROUTINE falfluxhucm_z

   SUBROUTINE onecond1 &
                       & (tt,qq,pp,ror,vr1        &
                       & ,del1n,del2n,div1,div2   &
                       & ,ff1,psi1,r1,rlec,ro1bl  &
                       & ,aa1_my,bb1_my,aa2_my,bb2_my &
                       & ,col,dtcond,icemax,nkr,isym1 &
                       & ,isym2,isym3,isym4,isym5,iin,kin,w_in,ccn_reg &
                       & ,latheatfac1,latheatfac2,latheatfac3,latheatfac4,supsat_out,sat_method,rho)
     IMPLICIT NONE
     REAL(KIND=wp),INTENT(INOUT) :: supsat_out
     REAL(KIND=wp),INTENT(IN) :: rho
     INTEGER, INTENT(IN) :: sat_method
     LOGICAL,  INTENT(IN)  :: &
      & latheatfac1, latheatfac2, latheatfac3, latheatfac4
     INTEGER nkr,icemax, isym1, isym2(icemax),isym3,isym4,isym5, iin, kin
     REAL(KIND=wp) ::   col,vr1(nkr)         &
      &       ,aa1_my,bb1_my,aa2_my,bb2_my &
      &       ,dtcond, w_in,ccn_reg 
     INTEGER i_bergeron, & !i_abergeron
      & kr,itime,kcond,nr,nrm
     REAL(KIND=wp) :: al1,al2,d,gam,pod, &
      & rv_my,d_myin, & !dtlref, dt0lref
!     & a1_myn, bb1_myn, a2_myn, bb2_myn, &
      & dt,dtt,xrad, &
      & epsdel, epsdel2,&
      & ror, &
      & cwhucm,b6,b8l, & !b8i, &
      & del1,del2,del1s,del2s, &
      & timenew,timerev,sfn11,sfn12, &
      & sfnl,sfni,b5l,b5i,b7l,b7i,dopl,dopi,rw,ri,qw,pw, &
      & pi,dtnewl,d1n,d2n
     REAL(KIND=wp) :: dt_water_cond,dt_water_evap
     INTEGER k
     REAL(KIND=wp) :: ff1_old(nkr),supintw(nkr)
     DOUBLE PRECISION dsupintw(nkr),dd1n,db11_my,dal1
     DOUBLE PRECISION col3,rori,tpn,tps,qpn,qps,told,qold &
      &                  ,fi1_k &
      &                  ,r1_k  &
      &                  ,fi1r1 &
      &                  ,rmasslaa,rmasslbb     &
      &                  ,es1n,es2n,ew1n,argexp &
      &                  ,tt,qq,pp &
      &                  ,del1n,del2n,div1,div2 &
      &                  ,oper2,oper3,ar1,ar2
     DOUBLE PRECISION delmassl1

     ! droplets
     REAL(KIND=wp) :: r1(nkr) &
      &           ,rlec(nkr),ro1bl(nkr) &
      &           ,fi1(nkr),ff1(nkr),psi1(nkr) &
      &           ,b11_my(nkr) !,b12_my(nkr)

     ! new algorithm of mixed phase for evaporation

     ! new algorithm (no type of ice)
     REAL(KIND=wp) :: fl1(nkr),sfndummy(3),totccn_before, totccn_after
     INTEGER :: idrop
     DOUBLE PRECISION :: r1d(nkr),r1nd(nkr)
     oper2(ar1)=0.622/(0.622+0.378*ar1)/ar1
     oper3(ar1,ar2)=ar1*ar2/(0.622+0.378*ar1)
     DATA al1 /2500./, al2 /2834./, d /0.211/ &
      &      ,gam /1.e-4/, pod /10./
     DATA rv_my,d_myin &
      &      /461.5,0.211e-4/
!     DATA a1_myn, bb1_myn, a2_myn, bb2_myn &
!      &      /2.53,5.42,3.41e1,6.13/
     DATA epsdel, epsdel2 /0.1e-03,0.1e-03/
     DOUBLE PRECISION :: del1_d , del2_d, rw_d , pw_d, ri_d, pi_d, d1n_d, d2n_d, &
                            vr1_d(nkr)
     sfndummy = 0.0
!    b12_my = 0.0
     b11_my = 0.0
!    i_abergeron=0
     i_bergeron=0
     col3=3.0*col
     itime=0
     kcond=0
     dt_water_cond=0.4
     dt_water_evap=0.4
     itime=0
     kcond=0
!    dt0lref=0.2
!    dtlref=0.4

     nr=nkr
     nrm=nkr-1
     dt=dtcond
     dtt=dtcond
     xrad=0.
     cwhucm=0.
     xrad=0.
     b6=cwhucm*gam-xrad
     b8l=1./ror
!    b8i=1./ror
     rori=1./ror

     ! ... ccn_regeneration
     totccn_before = 0.0
     totccn_before = sum(psi1(1:nkr)*r1(1:nkr))*3.0*col
     DO kr=1,nkr
       ff1_old(kr)=ff1(kr)
       supintw(kr)=0.0
       dsupintw(kr)=0.0
     END DO

     tpn=tt
     qpn=qq
     DO kr=1,nkr
       fi1(kr)=ff1(kr)
     END DO

     ! warm mp (condensation or evaporation) (begin)
     timenew=0.
     itime=0
     told = tpn
     qold = qpn
     r1d = r1
     r1nd = r1d
     sfnl = 0.0
     sfn11 = 0.0

 56  itime = itime+1
     timerev = dt-timenew
     timerev = dt-timenew
     del1 = del1n
     del2 = del2n
     del1s = del1n
     del2s = del2n
     tps = tpn
     qps = qpn

     IF (isym1 == 1) THEN !water droplets exist
       fl1 = 0.0
       vr1_d = vr1
       ! calculation of f_1 in (3.10) khain&sednev,1996:
       CALL jerrate_ks(r1d,tps,pp,vr1_d,rlec,ro1bl,b11_my,1,1,fl1,nkr,icemax,latheatfac4,sat_method)
       sfndummy(1)=sfn11
       ! calculation of coefficients in equations for supersaturation
       CALL jertimesc_ks(fi1,r1d,sfndummy,b11_my,b8l,1,nkr,col)
       sfn11 = sfndummy(1)
     END IF

     sfn12 = 0.0
     sfnl = sfn11 + sfn12
     sfni = 0.
     b5l=bb1_my/tps/tps
     b5i=bb2_my/tps/tps
     b7l=b5l*b6
     b7i=b5i*b6
     dopl=1.+del1s
     dopi=1.+del2s

     ! after each substep, the new sup sat over water and ice are calculated. These values are used to 
     ! calculate the new temperature and humidity at each substep. Therefore sup sat over ice is calculated even
     ! when there is no ice. The procedure can in future be simplified by assuming a linear changes of T and Q

     IF (.not. latheatfac1) THEN
      rw=(oper2(qps)+b5l*al1)*dopl*sfnl !r1 coeff in (3.11) in khain&sednev,1996, al1=lw/cp
     ELSE
      rw=(oper2(qps)+b5l*al1*cpd/cvd)*dopl*sfnl !r1 coeff in (3.11) in khain&sednev,1996, al1=lw/cp
     END IF
     ri=(oper2(qps)+b5l*al2)*dopl*sfni !r2 coeff in (3.11) in khain&sednev,1996, al2=li/cp
     qw=b7l*dopl
     IF (.not. latheatfac1) THEN
       pw=(oper2(qps)+b5i*al1)*dopi*sfnl !p1 coeff in (3.11) in khain&sednev,1996
     ELSE
       pw=(oper2(qps)+b5i*al1*cpd/cvd)*dopi*sfnl !p1 coeff in (3.11) in khain&sednev,1996
     END IF
     pi=(oper2(qps)+b5i*al2)*dopi*sfni !p2 coeff in (3.11) in khain&sednev,1996
!    qi=b7i*dopi
     IF (rw.NE.rw .OR. pw.NE.pw)THEN
       PRINT*, 'nan in onecond1'
       CALL finish(TRIM(modname),"fatal error in onecond1 (rw or pw are nan), model stop")
     END IF

     ! every ncond substep can be still too large for droplets growth. we want that during diffusional growth, 
     ! the droplet radius change will be less then few bins. therefore ncond substep is further divided by up 
     ! to kcond=10. in case that after 10 steps, time will be still lower than ncond sub step, the last step 
     ! will have this delta time. note that this is done without updates of t,q,s
     kcond=10 ! kcond is a flag
     IF (del1n >= 0.0d0) kcond=11 ! is it diffusional growth or evaporation??? --> kcond
     IF (kcond == 11) THEN
       dtnewl = dt
       dtnewl = MIN(dtnewl,timerev) ! timerev are small substeps within ncond sub step
       timenew = timenew + dtnewl
       dtt = dtnewl
       IF (dtt < 0.0) CALL finish(TRIM(modname),"fatal error in onecond1-del1n>0:(dtt<0), model stop")
       del1_d = del1
       del2_d = del2
       rw_d = rw
       pw_d = pw
       ri_d = ri
       pi_d = pi

       ! solving the equation for 2 supersaturations del1n,del2n
       CALL jersupsat_ks(del1_d,del2_d,del1n,del2n, &
                  rw_d,pw_d,ri_d,pi_d, &
                  dtt,d1n_d,d2n_d,0.0_wp,0.0_wp, &
                  isym1,isym2,isym3,isym4,isym5)
       del1 = del1_d
       del2 = del2_d
       rw = rw_d
       pw = pw_d
       ri = ri_d
       pi = pi_d
       d1n = d1n_d
       d2n = d2n_d

       IF (isym1 == 1)THEN ! water exists
         idrop = isym1
         !bin mass change and remapping - (3.14) in khain&sednev,1996:
         CALL jerdfun_ks(r1d, r1nd, b11_my, fi1, psi1, d1n, &
                    isym1, 1, 1, idrop, nkr, col, 1, iin, kin)
       END IF

       IF ((del1.GT.0.and.del1n.LT.0).and.abs(del1n).GT.epsdel) THEN
         CALL finish(TRIM(modname),"fatal error in onecond1-1 (del1.GT.0.and.del1n.LT.0), model stop")
       END IF
     ELSE
       ! evaporation - only water
       ! in case : kcond.NE.11
!      dtimeo = dt
       dtnewl = dt
       dtnewl = MIN(dtnewl,timerev)
       timenew = timenew + dtnewl
       dtt = dtnewl
       IF (dtt < 0.0) CALL finish(TRIM(modname),"fatal error in onecond1-del1n<0:(dtt<0), model stop")
       del1_d = del1
       del2_d = del2
       rw_d = rw
       pw_d = pw
       ri_d = ri
       pi_d = pi

       ! solving the equation for 2 supersaturations del1n,del2n. Check, why division of the code for diff 
       ! growth and evap is needed (???):
       CALL jersupsat_ks(del1_d,del2_d,del1n,del2n, &
                 rw_d,pw_d,ri_d,pi_d, &
                 dtt,d1n_d,d2n_d,0.0_wp,0.0_wp, &
                 isym1,isym2,isym3,isym4,isym5)
       del1 = del1_d
       del2 = del2_d
       rw = rw_d
       pw = pw_d
       ri = ri_d
       pi = pi_d
       d1n = d1n_d
       d2n = d2n_d
       IF (isym1 == 1)THEN
         idrop = isym1
         !bin mass change and remapping - (3.14) in khain&sednev,1996:
         CALL jerdfun_ks(r1d, r1nd, b11_my, &
         fi1, psi1, d1n, &
         isym1, 1, 1, idrop, nkr, col, 1, iin, kin)
       END IF
       IF ((del1.LT.0.and.del1n.GT.0).and.abs(del1n).GT.epsdel) THEN
         CALL finish(TRIM(modname),"fatal error in onecond1-2 (del1.LT.0.and.del1n.GT.0), model stop")
       END IF
     END IF

     ! masses:
     rmasslbb=0.
     rmasslaa=0.

     ! ... before jernewf (only water)
     DO k=1,nkr
       fi1_k = fi1(k)
       r1_k = r1(k) !mass of bin k, named xl outside the SUBROUTINE
       fi1r1 = fi1_k*r1_k*r1_k
       rmasslbb = rmasslbb+fi1r1 !calculate cloud water content before diffusional growth ncond substep
     END DO
     rmasslbb = rmasslbb*col3*rori !result: cloud water content before diffusional growth ncond substep
     IF (rmasslbb.LE.0.) rmasslbb=0.
     ! ... after jernewf (only water)
     DO k=1,nkr
       fi1_k=psi1(k)
       r1_k=r1(k)
       fi1r1=fi1_k*r1_k*r1_k
       rmasslaa=rmasslaa+fi1r1
     END DO
     rmasslaa=rmasslaa*col3*rori !result: cloud water content after diffusional growth ncond substep
     IF (rmasslaa.LE.0.) rmasslaa=0.

     delmassl1 = rmasslaa - rmasslbb !cloud mass change during substep (new-old)
     qpn = qps - delmassl1 !new specific humidity
     dal1 = al1 !L/cp
     IF (.not. latheatfac2) THEN
       tpn = tps + dal1*delmassl1 ! new temperature after substep (latent heat release or evaporation)
     ELSE
       tpn = tps + dal1*delmassl1*cpd/cvd ! new temperature after substep (latent heat release or evaporation)
     END IF

     IF (abs(dal1*delmassl1) > 10.0 )THEN !temperature change during one substep (latent heat release or evaporation)
       PRINT*,"onecond1-in(start)"
       PRINT*,"i=",iin,"kin",kin,"w",w_in 
       PRINT*,"delmassl1",delmassl1,"dt",dtt
       PRINT*,"del1n,del2n,del1,del2,d1n,d2n,rw,pw,ri,pi,dt"
       PRINT*,del1n,del2n,del1,del2,d1n,d2n,rw,pw,ri,pi,dtt
       PRINT*,"tps",tps,"qps",qps
       PRINT*,'fi1 before',fi1,'psi1 after',psi1
       PRINT*,"onecond1-in(end)"
       CALL finish(TRIM(modname),"fatal error in onecond1-in (abs(dal1*delmassl1) > 3.0), model stop")
     END IF

     ! ... supersaturation (only water)
     IF (sat_method == 1) THEN
       argexp=-bb1_my/tpn
       es1n=aa1_my*EXP(argexp)
       argexp=-bb2_my/tpn
       es2n=aa2_my*EXP(argexp)
     ELSE
       es1n=10.0*610.78*EXP(17.269*(tpn-273.15)/(tpn-35.86))
       es2n=10.0*610.78*EXP(21.875*(tpn-273.15)/(tpn-7.66))
     END IF

     IF ((sat_method == 1) .OR. (sat_method == 2)) THEN
       ew1n=oper3(qpn,pp)
     ELSE
       ew1n=10.0*qpn*rho*r_v*tpn !icon
     END IF

     IF (es1n == 0.0d0) THEN
       del1n=0.5
       div1=1.5
       !print*,'es1n onecond1 = 0'
       !call finish(trim(modname),"fatal error in onecond1 (es1n.eq.0), model stop")
     ELSE
       div1 = ew1n/es1n
       del1n = ew1n/es1n-1.
     END IF
     IF (es2n == 0.0d0)THEN
       del2n=0.5
       div2=1.5
       !print*,'es2n onecond1 = 0'
       !call finish(trim(modname),"fatal error in onecond1 (es2n.eq.0), model stop")
     ELSE
       del2n = ew1n/es2n-1.
       div2 = ew1n/es2n
     END IF
     ! calculation of the full integral over the sup sat (dsupintw), without ncond substeps, and 
     ! then perform a single remapping, instead of doing remapping after each substep
     IF (isym1 == 1) THEN
       DO kr=1,nkr
         supintw(kr)=supintw(kr)+b11_my(kr)*d1n
         dd1n=d1n
         db11_my=b11_my(kr)
         dsupintw(kr)=dsupintw(kr)+db11_my*dd1n
       END DO
     END IF

     ! ... repeate time step (only water: condensation or evaporation)
     IF (timenew.LT.dt) goto 56

!!57   CONTINUE

     IF (isym1 == 1) THEN
       ! a single diffusional growth step using the integral of sup sat. we do it for water only, 
       ! sinceexact collision initiation time is important, to prevent water psd broadening and early rain initiation.
       !bin mass change and remapping - (3.14) in khain&sednev,1996:
       !the difference between jerdfun_ks and jerdfun_new_ks is: 
       !in jerdfun_ks (3.14) is summed up within the ncond loop, and the remapping is done ncond times, but only for t update and not for the update of the final psd.
       !in jerdfun_new_ks (only for water) (3.14) is summed up at the end, after ncond loop, so that the remapping which updates the final psd is performed only once:
       CALL jerdfun_new_ks (r1d,r1nd,supintw, &
             ff1_old,psi1, &
             idrop, nkr, col,1,iin,kin)
     END IF

     rmasslaa=0.0
     rmasslbb=0.0

     DO k=1,nkr
       fi1_k=ff1_old(k)
       r1_k=r1(k)
       fi1r1=fi1_k*r1_k*r1_k
       rmasslbb=rmasslbb+fi1r1
     END DO
     rmasslbb=rmasslbb*col3*rori
     IF (rmasslbb.LT.0.0) rmasslbb=0.0

     DO k=1,nkr
       fi1_k=psi1(k)
       r1_k=r1(k)
       fi1r1=fi1_k*r1_k*r1_k
       rmasslaa=rmasslaa+fi1r1
     END DO
     rmasslaa=rmasslaa*col3*rori
     IF (rmasslaa.LT.0.0) rmasslaa=0.0
     delmassl1 = rmasslaa-rmasslbb

     !latent heat release:
     qpn = qold - delmassl1
     dal1 = al1
     IF (.not. latheatfac3) THEN
       tpn = told + dal1*delmassl1
       lh_ce_1 = lh_ce_1 + dal1*delmassl1
     ELSE
       tpn = told + dal1*delmassl1*cpd/cvd
       lh_ce_1 = lh_ce_1 + dal1*delmassl1*cpd/cvd
     END IF

     ! ... ccn regeneration
     totccn_after = 0.0_wp
     totccn_after = sum(psi1(1:nkr)*r1(1:nkr))*3.0*col ! [cm-3]
     ccn_reg = ccn_reg + max((totccn_before - totccn_after),0.0_wp)
     IF (abs(dal1*delmassl1) > 10.0_wp )THEN
       PRINT*,"onecond1-out (start)"
       PRINT*,"i=",iin,"kin",kin,"w",w_in 
       PRINT*,"del1n,del2n,d1n,d2n,rw,pw,ri,pi,dt"
       PRINT*,del1n,del2n,d1n,d2n,rw,pw,ri,pi,dtt
       PRINT*,"i=",iin,"kin",kin
       PRINT*,"tps=",tps,"qps=",qps,"delmassl1",delmassl1
       PRINT*,"dal1=",dal1
       PRINT*,rmasslbb,rmasslaa
       PRINT*,"fi1",fi1
       PRINT*,"psi1",psi1
       PRINT*,"onecond1-out (end)"
       IF (abs(dal1*delmassl1) > 5.0 )THEN
         CALL finish(TRIM(modname),"fatal error in onecond1-out (abs(dal1*delmassl1) > 5.0), model stop")
       END IF
     END IF

     ! ... supersaturation
     IF (sat_method == 1) THEN
       argexp=-bb1_my/tpn
       es1n=aa1_my*EXP(argexp)
       argexp=-bb2_my/tpn
       es2n=aa2_my*EXP(argexp)
     ELSE
       es1n=10.0*610.78*EXP(17.269*(tpn-273.15)/(tpn-35.86))
       es2n=10.0*610.78*EXP(21.875*(tpn-273.15)/(tpn-7.66))
     END IF

     IF ((sat_method == 1) .OR. (sat_method == 2)) THEN
       ew1n=oper3(qpn,pp)
     ELSE
       ew1n=10.0*qpn*rho*r_v*tpn !icon
     END IF

     IF (es1n == 0.0d0) THEN
       del1n=0.5
       div1=1.5
       CALL finish(TRIM(modname),"fatal error in onecond1 (es1n.EQ.0), model stop")
     ELSE
       div1=ew1n/es1n
       del1n=ew1n/es1n-1.
     END IF
     IF (es2n.EQ.0)THEN
       del2n=0.5
       div2=1.5
       CALL finish(TRIM(modname),"fatal error in onecond1 (es2n.EQ.0), model stop")
     ELSE
       del2n=ew1n/es2n-1.
       div2=ew1n/es2n
     END IF

     tt=tpn
     qq=qpn
     DO kr=1,nkr
       ff1(kr)=psi1(kr)
     END DO
     supsat_out=del1n 
     
     RETURN
   END SUBROUTINE onecond1

   SUBROUTINE coal_bott_new_warm(ff1r, tt, pp, dt_coll, krdrop)                  
     IMPLICIT NONE

     INTEGER,INTENT(IN) :: krdrop
     REAL(KIND=wp),INTENT(IN) :: dt_coll
     REAL(KIND=wp),INTENT(INOUT) :: ff1r(:)
     REAL(KIND=wp),INTENT(INOUT) :: tt
     REAL(KIND=wp),INTENT(IN) :: pp
     INTEGER :: kr,icol_drop,icol_drop_brk
     REAL(KIND=wp) :: g1(nkr), &
                      gdumb(jmax),gdumb_bf_breakup(jmax),xl_dumb(jmax)
     REAL(KIND=wp) :: t_new, pp_r
     INTEGER :: it, ndiv, it_is_rain1, it_is_cloud
     REAL(KIND=wp) :: break_drop_bef,break_drop_aft,dtbreakup,break_drop_per,   &
                      fl1(nkr),                                        &
                      cld_dsd_nbf,cld_dsd_naf,cld_dsd_mbf, & 
                      cld_dsd_maf,rain_dsd_nbf, &
                      rain_dsd_naf,rain_dsd_maf,rain_dsd_mbf,dsd_hlp(nkr), &
                      mass_flux(nkr)
     REAL(KIND=wp), PARAMETER :: g_lim = 1.0d-19*1.0d3

     icol_drop_brk=0
     icol_drop=0

     t_new = tt
     pp_r = pp

! probability of collision between 2 particles per unit of time, including swept volume and collision efficiency:
! Note that collision kernals were calculated at 1000, 750 and 500 mb and here it will be intepolated to height levels:
     CALL kernals_ks(dt_coll,nkr,pp_r)

! in case of ice, sticking efficiency depends on temperature
! and modkrn_ks is called. Here - not needed

! the calculation of collisions will be done for distribution funcions with respect to masses (mg)
! and therefore we do here translation from PSD (integral=concentration) to PSD (integral==mass content):

     DO kr=1,nkr !loop over all bins
       g1(kr) = ff1r(kr)*3.0*xl(kr)*xl(kr)*1.0e3        ! *1.e3 for g->mg  

       ! check whether collision/breakup is possible for each type, i.e. 
       ! if at least one bin is non-zero, we will call collision for this type:
       IF (kr > krmin_breakup .and. g1(kr) > g_lim) icol_drop_brk = 1
       IF (ibreakup == 0) icol_drop_brk = 0
       IF (g1(kr) > g_lim) icol_drop=1
     END DO

! ----------------------------
! ... drop-drop collisions
! ----------------------------

     !set icol_drop=0 to close collisions
     IF (icol_drop == 1) THEN

       ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       ! here we calculate the autoconversion rate using a dummy call to coll_xxx_lwf();     
       ! evaluating autoconv. from the cloud-mode spectra
    
       dsd_hlp = g1
       rain_dsd_nbf = col*sum(dsd_hlp(krdrop+1:nkr)/xl(krdrop+1:nkr))*1.e-3                                          ! in [#/cm3]
       rain_dsd_mbf = col*sum(dsd_hlp(krdrop+1:nkr))*1.e-3                                                           ! in [g/cm3]
       cld_dsd_mbf  = col*sum(dsd_hlp(1:krdrop))*1.e-3                                                               ! in [g/cm3]
       cld_dsd_nbf  = col*sum(dsd_hlp(1:krdrop)/xl(1:krdrop))*1.e-3                                                  ! in [#/cm3]
       mass_flux = 0.0                                                  
       ! x1 - the larger particle type, x2 - the smaller particle type, x3 - the resulting particle type
       ! therefore "xxx" means water + water --> water, or snow + snow --> snow

       CALL coll_xxx_bott_mod1(dsd_hlp,1,krdrop,krdrop,cwll,xl_mg,chucm,ima,1.0d0,nkr,mass_flux)

       rain_dsd_naf = col*sum(dsd_hlp(krdrop+1:nkr)/xl(krdrop+1:nkr))*1.e-3                                          ! in [#/cm3]
       rain_dsd_maf = col*sum(dsd_hlp(krdrop+1:nkr))*1.e-3                                                           ! in [g/cm3]
       cld_dsd_maf  = col*sum(dsd_hlp(1:krdrop))*1.e-3                                                               ! in [g/cm3]
       cld_dsd_naf  = col*sum(dsd_hlp(1:krdrop)/xl(1:krdrop))*1.e-3                                                  ! in [#/cm3]
       it_is_cloud = 0
       IF (cld_dsd_mbf > 0.01*1.0e-6) it_is_cloud = 1 
       IF ( it_is_cloud == 1 ) THEN
         auto_cld_msink_b  = rain_dsd_maf - rain_dsd_mbf                                                              ! [+]
         auto_cld_nsink_b  = cld_dsd_nbf  - cld_dsd_naf                                                               ! [+]  
       END IF

       ! evaluating accretion from the full spectra
       dsd_hlp = g1
       cld_dsd_nbf  = col*sum(dsd_hlp(1:krdrop)/xl(1:krdrop))*1.e-3                                                   ! in [#/cm3]
       cld_dsd_mbf  = col*sum(dsd_hlp(1:krdrop))*1.e-3                                                                ! in [g/cm3]
       rain_dsd_nbf = col*sum(dsd_hlp(krdrop+1:nkr)/xl(krdrop+1:nkr))*1.e-3                                           ! in [#/cm3]
       rain_dsd_mbf = col*sum(dsd_hlp(krdrop+1:nkr))*1.e-3                                                            ! in [g/cm3]
       mass_flux = 0.0

       CALL coll_xxx_bott_mod2(dsd_hlp,1,krdrop,krdrop+1,nkr,cwll,xl_mg,chucm,ima,1.0d0,nkr,mass_flux)

       cld_dsd_naf  = col*sum(dsd_hlp(1:krdrop)/xl(1:krdrop))*1.e-3                                                   ! in [#/cm3]
       cld_dsd_maf  = col*sum(dsd_hlp(1:krdrop))*1.e-3                                                                ! in [g/cm3]
       rain_dsd_naf = col*sum(dsd_hlp(krdrop+1:nkr)/xl(krdrop+1:nkr))*1.e-3                                           ! in [#/cm3]
       rain_dsd_maf = col*sum(dsd_hlp(krdrop+1:nkr))*1.e-3                                                            ! in [g/cm3]
       it_is_cloud = 0
       IF (cld_dsd_mbf > 0.01*1.0e-6) it_is_cloud = 1
       IF ( it_is_cloud == 1 ) THEN
         accr_cld_nsink_b  = cld_dsd_nbf - cld_dsd_naf                                                                ! [+]   
         accr_cld_msink_b  = cld_dsd_mbf - cld_dsd_maf                                                                ! [+]
       END IF

       ! evaluating rain self-collection from the rain spectra
       dsd_hlp = g1
       rain_dsd_mbf = col*sum(dsd_hlp(krdrop+1:nkr))*1.e-3
       rain_dsd_nbf = col*sum(dsd_hlp(krdrop+1:nkr)/xl(krdrop+1:nkr))*1.e-3
       mass_flux = 0.0                                         
       it_is_rain1 = 0 

       CALL coll_xxx_bott_mod1(dsd_hlp,krdrop+1,nkr,nkr,cwll,xl_mg,chucm,ima,1.0d0,nkr,mass_flux)

       rain_dsd_naf = col*sum(dsd_hlp(krdrop+1:nkr)/xl(krdrop+1:nkr))*1.e-3                                           ! in [#/cm3]
       IF (rain_dsd_mbf > 1.0d-50) it_is_rain1 = 1
       IF ( it_is_rain1 == 1 ) THEN
         selfc_rain_nchng_b = max(0.0d0,rain_dsd_nbf - rain_dsd_naf)
       END IF
       ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

       ! here we call the standard 'coll_xxx()' routine with the full spectrum (dsd+rsd)                              ! in [g/cm3]
       fl1 = 1.0
       mass_flux = 0.d0
       ! input:  g1 - mass PSD, cwll - collision kernel liquid-liquid, xl_mg - mass bins:
       ! output: new g1:
       CALL coll_xxx_bott(g1,cwll,xl_mg,chucm,ima,1.0d0,nkr,mass_flux)                                        

! --------------------------------------------------------
! ... probability of drop breakup after collision:
! --------------------------------------------------------

! There is an impirical distribution of fragments sizes after breakup
! Here we do iterative process to use this distribution with the restriction of mass conservation:
       IF (icol_drop_brk == 1)THEN
         ndiv = 1
10       CONTINUE
         DO it = 1,ndiv             ! ndiv=number of iterations
           dtbreakup = dt_coll/ndiv ! timestep of brekup
           IF (it == 1)THEN
             DO kr=1,jmax
               gdumb(kr)= g1(kr)*1.d-3
               gdumb_bf_breakup(kr) =  g1(kr)*1.d-3
               xl_dumb(kr)=xl_mg(kr)*1.d-3
             END DO
             break_drop_bef=0.d0
             DO kr=1,jmax
               break_drop_bef = break_drop_bef+g1(kr)*1.d-3
             END DO
           END IF
           CALL coll_breakup_ks(gdumb, xl_dumb, jmax, dtbreakup, jbreak, pkij, qkj, nkr, nkr)
         END DO

! check for errors:
         DO kr=1,nkr
           ff1r(kr) = (1.0d3*gdumb(kr))/(3.0*xl(kr)*xl(kr)*1.e3)
           IF (ff1r(kr) < 0.0)THEN
             IF (ndiv < 8)THEN
               ndiv = 2*ndiv   ! in this case calculate breakup again with twice more iterations
               go to 10
             ELSE              ! if ndiv is already large and still the final PSD is <0, give up:
               !PRINT*,"nobreakup",jin,kin,itimestep,ndiv
               go to 11
             END IF
           END IF
           IF (ff1r(kr) .NE. ff1r(kr)) THEN ! check for nan
             PRINT*,kr,gdumb(kr),gdumb_bf_breakup(kr),xl(kr)
             PRINT*,it,ndiv, dtbreakup
             PRINT*,gdumb
             PRINT*,gdumb_bf_breakup
             CALL finish(TRIM(modname),"in coal_bott after coll_breakup - ff1r nan, model stop")
           END IF
         END DO

! mass conservation:
         break_drop_aft=0.0d0
         DO kr=1,jmax
           break_drop_aft=break_drop_aft+gdumb(kr)
         END DO
         break_drop_per=break_drop_aft/break_drop_bef
         IF (break_drop_per > 1.001)THEN
           ndiv=ndiv*2  ! in this case calculate breakup again with twice more iterations
           go to 10
         ELSE
           DO kr=1,jmax
             g1(kr) = gdumb(kr)*1.d3
           END DO
         END IF
       END IF ! IF icol_drop_brk.EQ.1
     END IF ! IF icol_drop.EQ.1

11   CONTINUE

     ! recalculation of density FUNCTION f1,f3,f4,f5 in  units [1/(g*cm**3)] :
     DO kr=1,nkr
       ff1r(kr)=g1(kr)/(3.*xl(kr)*xl(kr)*1.e3)
       IF ((ff1r(kr) .NE. ff1r(kr)) .OR. ff1r(kr) < 0.0)THEN
         PRINT*,"g1",g1
         CALL finish(TRIM(modname),"stop at end coal_bott - ff1r nan or ff1r < 0.0, model stop")
       END IF
     END DO

     IF (abs(tt-t_new).GT.5.0) THEN
       CALL finish(TRIM(modname),"fatal error in module_mp_warm_sbm del_t 5 k, model stop")
     END IF
  
     tt = t_new

     RETURN
   END SUBROUTINE coal_bott_new_warm

   SUBROUTINE coll_xxx_lwf(g,fl,ckxx,x,c,ima,prdkrn,nkr,output_flux)
     IMPLICIT NONE
     INTEGER,INTENT(IN) :: nkr
     REAL(KIND=wp),INTENT(INOUT) :: g(:),fl(:)
     REAL(KIND=wp),INTENT(IN) :: ckxx(:,:),x(:), c(:,:)
     INTEGER,INTENT(IN) :: ima(:,:)
     REAL(KIND=wp),INTENT(IN) :: prdkrn
     REAL(KIND=wp),INTENT(INOUT) :: output_flux(:)
     REAL(KIND=wp):: gmin,x01,x02,x03,gsi,gsj,gsk,gsi_w,gsj_w,gsk_w,gk, &
                          gk_w,fl_gk,fl_gsk,flux,x1,flux_w,g_k_w,g_kp_old,g_kp_w
     INTEGER :: i,ix0,ix1,j,k,kp
     
     gmin=g_lim*1.0d3
   
   ! ix0 - lower limit of integration by i
     DO i=1,nkr-1
       ix0=i
       IF (g(i).GT.gmin) goto 2000
     END DO
2000 CONTINUE
     IF (ix0.EQ.nkr-1) RETURN
   
   ! ix1 - upper limit of integration by i
     DO i=nkr-1,1,-1
       ix1=i
       IF (g(i).GT.gmin) goto 2010
     END DO
2010 CONTINUE

   ! ... collisions
     DO i=ix0,ix1
       IF (g(i).LE.gmin) goto 2020
       DO j=i,ix1
         IF (g(j).LE.gmin) goto 2021
         k=ima(i,j)
         kp=k+1
         x01=ckxx(i,j)*g(i)*g(j)*prdkrn
         x02=dmin1(x01,g(i)*x(j))
         IF (j.NE.k) x03=dmin1(x02,g(j)*x(i))
         IF (j.EQ.k) x03=x02
         gsi=x03/x(j)
         gsj=x03/x(i)
         gsk=gsi+gsj
         IF (gsk.LE.gmin) goto 2021
         gsi_w=gsi*fl(i)
         gsj_w=gsj*fl(j)
         gsk_w=gsi_w+gsj_w
         gsk_w=dmin1(gsk_w,gsk)
         g(i)=g(i)-gsi
         g(i)=dmax1(g(i),0.0d0)
         g(j)=g(j)-gsj
         IF (j.NE.k) g(j)=dmax1(g(j),0.0d0)
         gk=g(k)+gsk
   
         IF (g(j).LT.0.d0.and.gk.LE.gmin) THEN
           g(j)=0.d0
           g(k)=g(k)+gsi
           goto 2021
         END IF
   
         IF (gk.LE.gmin) goto 2021
         gk_w=g(k)*fl(k)+gsk_w
         gk_w=dmin1(gk_w,gk)
         fl_gk=gk_w/gk
         fl_gsk=gsk_w/gsk
         flux=0.d0
         x1=DLOG(g(kp)/gk+1.d-15)
         flux=gsk/x1*(EXP(0.5d0*x1)-EXP(x1*(0.5d0-c(i,j))))
         flux=dmin1(flux,gsk)
         flux=dmin1(flux,gk)
         ! changed to >= and corrected to bin #33               
         IF (kp.GE.kp_flux_max) flux=0.5d0*flux
         flux_w=flux*fl_gsk
         flux_w=dmin1(flux_w,gsk_w)
         flux_w=dmin1(flux_w,gk_w)
         g(k)=gk-flux
         g(k)=dmax1(g(k),gmin)
         g_k_w=gk_w-flux_w
         g_k_w=dmin1(g_k_w,g(k))
         g_k_w=dmax1(g_k_w,0.0d0)
         fl(k)=g_k_w/g(k)
         g_kp_old=g(kp)
         g(kp)=g(kp)+flux
         ! output flux - for autoconv.
         output_flux(kp) = output_flux(kp) + flux
         g(kp)=dmax1(g(kp),gmin)
         g_kp_w=g_kp_old*fl(kp)+flux_w
         g_kp_w=dmin1(g_kp_w,g(kp))
         fl(kp)=g_kp_w/g(kp)

         IF (fl(k).GT.1.0d0.and.fl(k).LE.1.0001d0) fl(k)=1.0d0
         IF (fl(kp).GT.1.0d0.and.fl(kp).LE.1.0001d0) fl(kp)=1.0d0
         IF (fl(k).GT.1.0001d0.OR.fl(kp).GT.1.0001d0 &
                 .OR.fl(k).LT.0.0d0.OR.fl(kp).LT.0.0d0) THEN
           PRINT*,    'in SUBROUTINE coll_xxx_lwf'
           PRINT*,    'snow - snow = snow'  
           IF (fl(k).GT.1.0001d0)  PRINT*, 'fl(k).GT.1.0001d0'
           IF (fl(kp).GT.1.0001d0) PRINT*, 'fl(kp).GT.1.0001d0'
           IF (fl(k).LT.0.0d0)  PRINT*, 'fl(k).LT.0.0d0'
           IF (fl(kp).LT.0.0d0) PRINT*, 'fl(kp).LT.0.0d0'
           PRINT*,    'i,j,k,kp'
           PRINT*,     i,j,k,kp
           PRINT*,    'ix0,ix1'
           PRINT*,     ix0,ix1   
           WRITE (*,'(A,4D13.5)') 'ckxx(i,j),x01,x02,x03:  ', ckxx(i,j),x01,x02,x03   
           WRITE (*,'(A,3D13.5)') 'gsi,gsj,gsk:  ', gsi,gsj,gsk   
           WRITE (*,'(A,3D13.5)') 'gsi_w,gsj_w,gsk_w:   ', gsi_w,gsj_w,gsk_w
           WRITE (*,'(A,2D13.5)') 'gk,gk_w:  ', gk,gk_w   
           WRITE (*,'(A,2D13.5)') 'fl_gk,fl_gsk:  ', fl_gk,fl_gsk   
           WRITE (*,'(A,2D13.5)') 'x1,c(i,j):  ', x1,c(i,j)   
           WRITE (*,'(A,D13.5)')  'flux:     ', flux
           WRITE (*,'(A,D13.5)')  'flux_w:   ', flux_w
           WRITE (*,'(A,D13.5)')  'g_k_w:    ', g_k_w
           WRITE (*,'(A,D13.5)')  'g_kP_w:   ', g_kp_w
           IF (fl(k).LT.0.0d0) PRINT*, &
                 'stop 2022: in SUBROUTINE coll_xxx_lwf, fl(k) < 0'
           IF (fl(kp).LT.0.0d0) PRINT*, &
                 'stop 2022: in SUBROUTINE coll_xxx_lwf, fl(kp) < 0'
           IF (fl(k).GT.1.0001d0) PRINT*, &
                 'stop 2022: in sub. coll_xxx_lwf, fl(k) > 1.0001'
           IF (fl(kp).GT.1.0001d0) PRINT*, &
                 'stop 2022: in sub. coll_xxx_lwf, fl(kp) > 1.0001'
              CALL finish(TRIM(modname),"in coal_bott sub. coll_xxx_lwf, model stop")
         END IF
2021     CONTINUE
       END DO
2020   CONTINUE
     END DO
   

     RETURN
   END SUBROUTINE coll_xxx_lwf

   SUBROUTINE coll_breakup_ks (gt_mg, xt_mg, jmax, dt, jbreak, &
                              pkij, qkj, nkrinput, nkr)
     IMPLICIT NONE
     INTEGER,INTENT(IN) :: jmax, jbreak, nkrinput, nkr
     REAL(KIND=wp),INTENT(IN) :: xt_mg(:), dt
     REAL(KIND=wp),INTENT(IN) :: pkij(:,:,:),qkj(:,:)
     REAL(KIND=wp),INTENT(INOUT) :: gt_mg(:)
     INTEGER :: ie, je, ke, nkrdiff, jdiff, k, i, j
     REAL(KIND=wp) :: gt(jmax), xt(jmax+1), ft(jmax), fa(jmax), dg(jmax), df(jmax), dbreak(jbreak) &
                    ,amweight(jbreak), gain, aloss
     ie=jbreak
     je=jbreak
     ke=jbreak
     ! gt_mg : mass distribution function of bott
     ! xt_mg : mass of bin in mg
     ! jmax  : number of bins
     ! dt    : timestep in s
     ! in cgs
     nkrdiff = nkrinput-nkr
     DO j=1,jmax
       xt(j)=xt_mg(j)
       gt(j)=gt_mg(j)
       ft(j)=gt(j)/xt(j)/xt(j)
     END DO
     !shift between coagulation and breakup grid
     jdiff=jmax-jbreak

     !initialization
     !shift to breakup grid
     fa = 0.0
     DO k=1,ke-nkrdiff
       fa(k)=ft(k+jdiff+nkrdiff)
     END DO

     !breakup: bleck's first order method
     !pkij: gain coefficients
     !qkj : loss coefficients

     xt(jmax+1)=xt(jmax)*2.0d0

     amweight = 0.0
     dbreak = 0.0
     DO k=1,ke-nkrdiff
       gain=0.0d0
       DO i=1,ie-nkrdiff
         DO j=1,i
           gain=gain+fa(i)*fa(j)*pkij(k,i,j)
         END DO
       END DO
       aloss=0.0d0
       DO j=1,je-nkrdiff
         aloss=aloss+fa(j)*qkj(k,j)
       END DO
       j=jmax-jbreak+k+nkrdiff
       amweight(k)=2.0/(xt(j+1)**2.0-xt(j)**2.0)
       dbreak(k)=amweight(k)*(gain-fa(k)*aloss)

       IF (dbreak(k) .NE. dbreak(k)) THEN
         PRINT*,dbreak(k),amweight(k),gain,fa(k),aloss
         PRINT*,"-"
         PRINT*,dbreak
         PRINT*,"-"
         PRINT*,amweight
         PRINT*,"-"
         PRINT*,j,jmax,jbreak,k,nkrdiff
         PRINT*,"-"
         PRINT*,fa
         PRINT*,"-"
         PRINT*,xt
         PRINT*,"-"
         PRINT*,gt
         CALL finish(TRIM(modname)," inside coll_breakup, nan, model stop")
       END IF
     END DO

     !shift rate to coagulation grid
     df = 0.0d0
     DO j=1,jdiff+nkrdiff
       df(j)=0.0d0
     END DO

     DO j=1,ke-nkrdiff
       df(j+jdiff)=dbreak(j)
     END DO

     !transformation to mass distribution function g(ln x)
     DO j=1,jmax
       dg(j)=df(j)*xt(j)*xt(j)
     END DO

     !time integration
     DO j=1,jmax
       gt(j)=gt(j)+dg(j)*dt
     END DO

     gt_mg = gt

     RETURN
   END SUBROUTINE coll_breakup_ks
  
   DOUBLE PRECISION FUNCTION polysvp (tt,itype,sat_method)

     IMPLICIT NONE
     REAL(KIND=wp), INTENT(IN) :: tt
     INTEGER,INTENT(IN) :: itype,sat_method
     REAL(wp), PARAMETER :: aa1_my = 2.53e12_wp, bb1_my = 5.42e3_wp, &
                            aa2_my = 3.41e13_wp, bb2_my = 6.13e3_wp
     REAL(wp) :: es1n, es2n
   
     method_select: select case(itype)
      ! liquid
      case(0)
        IF (sat_method == 1) THEN
          es1n = aa1_my*EXP(-bb1_my/tt)
        ELSE
          es1n=10.0*610.78*EXP(17.269*(tt-273.15)/(tt-35.86))
        END IF
        polysvp = es1n ! [dyn/cm2] to [mb]

      ! ice  
      case(1)
        IF (sat_method == 1) THEN
          es2n = aa2_my*EXP(-bb2_my/tt)
        ELSE
          es2n=10.0*610.78*EXP(21.875*(tt-273.15)/(tt-7.66))
        END IF
        polysvp = es2n ! [dyn/cm2] to [mb]

      case default
        polysvp = HUGE(1.0_dp)
    
      end select method_select

      RETURN
   END FUNCTION polysvp

   SUBROUTINE jerrate_ks (xls,tp,pp,vxl,riec,ro1bl,b11_my, &
                      id,in,fl1,nkr,icemax,latheatfac4,sat_method)
     IMPLICIT NONE
     LOGICAL,  INTENT(IN)  :: latheatfac4
     INTEGER,INTENT(IN) :: id, in, nkr, icemax,sat_method
     REAL(KIND=wp),INTENT(IN) :: ro1bl(nkr,id),riec(nkr,id),fl1(nkr)
     REAL(KIND=wp),INTENT(INOUT) :: b11_my(nkr,id)
     REAL(KIND=wp),INTENT(IN) :: pp, tp, xls(nkr,id),vxl(nkr,id)
     INTEGER :: kr, nskin(nkr), ice
     REAL(KIND=wp) :: ventplm(nkr), fd1(nkr,icemax),fk1(nkr,icemax), xl_my1(nkr,icemax), &
                          al1_my(2),esat1(2), tpreal
     REAL(KIND=wp) :: pzero, tzero, const, d_my, coeff_viscous, shmidt_number,     &
                          a, b, rvt, shmidt_number03, xls_kr_ice, ro1bl_kr_ice, vxl_kr_ice, reinolds_number, &
                          reshm, ventpl, constl, detl

     REAL(KIND=wp) :: deg01,deg03
     ! a1l_my - constants for "maxwell": mks
     REAL(KIND=wp), PARAMETER:: rv_my=461.5d4, cf_my=2.4d3, d_myin=0.211d0
     ! cgs :
     ! rv_my, cm*cm/sec/sec/kelvin - individual gas constant
     !                               for water vapour
     !rv_my=461.5d4
     ! d_myin, cm*cm/sec - coefficient of diffusion of water vapour
     !d_myin=0.211d0
     ! pzero, dynes/cm/cm - reference pressure
     pzero=1.013d6
     ! tzero, kelvin - reference temperature
     tzero=273.15d0

     DO kr=1,nkr
       IF (in==2 .and. fl1(kr)==0.0 .OR. in==6 .OR. in==3 .and. tp<273.15) THEN
         nskin(kr) = 2
       ELSE !in==1 or in==6 or lef/=0
         nskin(kr) = 1
       END IF
     END DO

     ! constants for clausius-clapeyron equation :
     ! a1_my(1),g/sec/sec/cm
     !	a1_my(1)=2.53d12
     ! a1_my(2),g/sec/sec/cm
     !	a1_my(2)=3.41d13
     ! bb1_my(1), kelvin
     !	bb1_my(1)=5.42d3
     ! bb1_my(2), kelvin
     !	bb1_my(2)=6.13d3
     ! al1_my(1), cm*cm/sec/sec - latent heat of vaporization
     al1_my(1)=2.5d10
     ! al1_my(2), cm*cm/sec/sec - latent heat of sublimation
     al1_my(2)=2.834d10
     ! cf_my, g*cm/sec/sec/sec/kelvin - coefficient of
     !                            thermal conductivity of air
     ! cf_my=2.4d3
     deg01=1.0/3.0
     deg03=1.0/3.0
     const=12.566372d0
     ! coefficient of diffusion
     d_my=d_myin*(pzero/pp)*(tp/tzero)**1.94d0
     ! coefficient of viscousity
     ! coeff_viscous=0.13 cm*cm/sec
     coeff_viscous=0.13d0
     ! shmidt number
     shmidt_number=coeff_viscous/d_my

     ! constants used for calculation of reinolds number

     a=2.0d0*(3.0d0/4.0d0/3.141593d0)**deg01
     b=a/coeff_viscous

     rvt=rv_my*tp
     ! update the saturation vapor pressure
     tpreal = tp
     esat1(1) = polysvp(tpreal,0,sat_method)
     esat1(2) = polysvp(tpreal,1,sat_method)
     DO kr=1,nkr
       ventplm(kr)=0.0d0
     END DO

     shmidt_number03=shmidt_number**deg03
     DO ice=1,id
       DO kr=1,nkr
         xls_kr_ice=xls(kr,ice)
         ro1bl_kr_ice=ro1bl(kr,ice)
         vxl_kr_ice=vxl(kr,ice)
         ! reynolds numbers
         reinolds_number= &
             b*vxl_kr_ice*(xls_kr_ice/ro1bl_kr_ice)**deg03
         reshm=dsqrt(reinolds_number)*shmidt_number03

         IF (reinolds_number<2.5d0) THEN
           ventpl=1.0d0+0.108d0*reshm*reshm
           ventplm(kr)=ventpl
         ELSE
           ventpl=0.78d0+0.308d0*reshm
           ventplm(kr)=ventpl
         END IF
       END DO

       ! ventpl_max is given in micro.prm include file
       DO kr=1,nkr
         ventpl=ventplm(kr)
         IF (ventpl>ventpl_max) THEN
           ventpl=ventpl_max
           ventplm(kr)=ventpl
         END IF
         constl=const*riec(kr,ice)
         fd1(kr,ice)=rvt/d_my/esat1(nskin(kr))
         IF (.not. latheatfac4) THEN
           fk1(kr,ice)=(al1_my(nskin(kr))/rvt-1.0d0)*al1_my(nskin(kr))/cf_my/tp
         ELSE
           fk1(kr,ice)=(al1_my(nskin(kr))*(cpd/cvd)/rvt-1.0d0)*al1_my(nskin(kr))*(cpd/cvd)/cf_my/tp
         END IF

         xl_my1(kr,ice)=ventpl*constl
         ! growth rate
         detl=fk1(kr,ice)+fd1(kr,ice)
         b11_my(kr,ice)=xl_my1(kr,ice)/detl
       END DO
     END DO

     RETURN
   END SUBROUTINE jerrate_ks

   SUBROUTINE jertimesc_ks (fi1,x1,sfn11, &
                            b11_my,cf,id,nkr,col)
     IMPLICIT NONE
     INTEGER,INTENT(IN) :: id,nkr
     REAL(KIND=wp),INTENT(IN) :: b11_my(nkr,id), fi1(nkr,id), col, cf
     REAL(KIND=wp),INTENT(IN) :: x1(nkr,id)
     REAL(KIND=wp),INTENT(OUT) :: sfn11(id)
     INTEGER :: ice, kr
     REAL(KIND=wp) :: sfn11s, fk, delm, fun, b11

     DO ice=1,id
       sfn11s=0.0d0
       sfn11(ice)=cf*sfn11s
       DO kr=1,nkr
         ! value of size distribution functions
         fk=fi1(kr,ice)
         ! delta-m
         delm=x1(kr,ice)*3.0d0*col
         ! integral's expression
         fun=fk*delm
         ! values of integrals
         b11=b11_my(kr,ice)
         sfn11s=sfn11s+fun*b11
       END DO
       ! cycle by kr
       ! correction
       sfn11(ice)=cf*sfn11s
     END DO
     ! cycle by ice
     RETURN
   END SUBROUTINE jertimesc_ks

   SUBROUTINE jersupsat_ks (del1,del2,del1n,del2n,         &
                            rw,pw,ri,pi,                   &
                            dt,del1int,del2int,dyn1,dyn2,  &
                            isym1,isym2,isym3,isym4,isym5)
     IMPLICIT NONE
     INTEGER,INTENT(INOUT) :: isym1, isym2(:), isym3, isym4, isym5
     REAL(KIND=wp),INTENT(IN) :: dt, dyn1, dyn2
     REAL(KIND=wp),INTENT(IN) :: del1, del2
     REAL(KIND=wp),INTENT(INOUT) :: del1n,del2n,del1int,del2int,rw, pw, ri, pi
     INTEGER ::  isymice, irw, iri !, ipi, ipw
     REAL(KIND=wp) :: x, expm1, deter, expr, expp, a, alfa, beta, gama, g31, g32, g2, expb, expg, &
                c11, c21, c12, c22, a1del1n, a2del1n, a3del1n, a4del1n, a1del1int, a2del1int, &
              a3del1int, a4del1int, a1del2n, a2del2n, a3del2n , a4del2n, a1del2int, a2del2int, &
              a3del2int, a4del2int, a5del2int

     expm1(x)=x+x*x/2.0d0+x*x*x/6.0d0+x*x*x*x/24.0d0+ &
                  x*x*x*x*x/120.0d0

     isymice = sum(isym2) + isym3 + isym4 + isym5
     irw = 1
!    ipw = 1
     iri = 1
!    ipi = 1  
     IF (MAX(rw,pw,ri,pi)<=rw_pw_ri_pi_min) THEN
       rw = 0.0
       irw = 0
       pw = 0.0
!      ipw = 0
       ri = 0.0
       iri = 0
       pi = 0.0
!      ipi = 0
       isym1 = 0
       isymice = 0
     ELSE
       IF (dmax1(rw,pw)>rw_pw_min) THEN
       ! a zero can pass through, assign a minimum value
         IF (rw < rw_pw_min*rw_pw_min) THEN
           rw = 1.0d-20
           irw = 0
         END IF
         IF (pw < rw_pw_min*rw_pw_min)THEN
           pw = 1.0d-20
!          ipw = 0
         END IF

         IF (dmax1(pi/pw,ri/rw)<=ratio_icew_min) THEN
           ! ... only water
           ri = 0.0
           iri = 0
           pi = 0.0
!          ipi = 0
           isymice = 0
         END IF

         IF (dmin1(pi/pw,ri/rw)>1.0d0/ratio_icew_min) THEN
           ! ... only ice
           rw = 0.0
           irw = 0
           pw = 0.0
!          ipw = 0
           isym1 = 0
         END IF
       ELSE
         ! only ice
         rw = 0.0
         irw = 0
         pw = 0.0
!        ipw = 0
         isym1 = 0
       END IF
     END IF

     IF (isymice == 0)THEN
       isym2 = 0
       isym3 = 0
       isym4 = 0
       isym5 = 0
     END IF

     deter=rw*pi-pw*ri
     IF (irw == 0 .and. iri == 0) THEN
       del1n=del1+dyn1*dt
       del2n=del2+dyn2*dt
       del1int=del1*dt+dyn1*dt*dt/2.0d0
       del2int=del2*dt+dyn2*dt*dt/2.0d0
goto   100
     END IF
     ! solution of equation for supersaturation with
     ! different deter values
     IF (iri == 0) THEN
       ! ... only water                                       (start)
       expr=EXP(-rw*dt)
       IF (abs(rw*dt)>1.0e-6) THEN
         del1n=del1*expr+(dyn1/rw)*(1.0d0-expr)
         del2n=pw*del1*expr/rw-pw*dyn1*dt/rw- &
               pw*dyn1*expr/(rw*rw)+dyn2*dt+ &
               del2-pw*del1/rw+pw*dyn1/(rw*rw)
         del1int=-del1*expr/rw+dyn1*dt/rw+ &
                 dyn1*expr/(rw*rw)+del1/rw-dyn1/(rw*rw)
         del2int=pw*del1*expr/(-rw*rw)-pw*dyn1*dt*dt/(2.0d0*rw)+ &
                 pw*dyn1*expr/(rw*rw*rw)+dyn2*dt*dt/2.0d0+ &
                 del2*dt-pw*del1*dt/rw+pw*dyn1*dt/(rw*rw)+ &
                 pw*del1/(rw*rw)-pw*dyn1/(rw*rw*rw)
goto     100
         ! in case dabs(rw*dt)>1.0d-6
       ELSE
         ! in case dabs(rw*dt)<=1.0d-6
         expr=expm1(-rw*dt)
         del1n=del1+del1*expr+(dyn1/rw)*(0.0d0-expr)
         del2n=pw*del1*expr/rw-pw*dyn1*dt/rw- &
               pw*dyn1*expr/(rw*rw)+dyn2*dt+del2
         del1int=-del1*expr/rw+dyn1*dt/rw+dyn1*expr/(rw*rw)
         del2int=pw*del1*expr/(-rw*rw)-pw*dyn1*dt*dt/(2.0d0*rw)+ &
               pw*dyn1*expr/(rw*rw*rw)+dyn2*dt*dt/2.0d0+ &
               del2*dt-pw*del1*dt/rw+pw*dyn1*dt/(rw*rw)
goto     100
       END IF
       ! ... only water                                                    (end)
       ! in case ri==0.0d0
     END IF

     IF (irw  ==  0) THEN
       ! ... only ice                                                    (start)
       expp=EXP(-pi*dt)
       IF (abs(pi*dt)>1.0e-6) THEN
         del2n = del2*expp+(dyn2/pi)*(1.0d0-expp)
         del2int = -del2*expp/pi+dyn2*dt/pi+ &
                  dyn2*expp/(pi*pi)+del2/pi-dyn2/(pi*pi)
         del1n = +ri*del2*expp/pi-ri*dyn2*dt/pi- &
                  ri*dyn2*expp/(pi*pi)+dyn1*dt+ &
                  del1-ri*del2/pi+ri*dyn2/(pi*pi)
         del1int = -ri*del2*expp/(pi*pi)-ri*dyn2*dt*dt/(2.0d0*pi)+ &
                  ri*dyn2*expp/(pi*pi*pi)+dyn1*dt*dt/2.0d0+ &
                  del1*dt-ri*del2*dt/pi+ri*dyn2*dt/(pi*pi)+ &
                  ri*del2/(pi*pi)-ri*dyn2/(pi*pi*pi)
goto     100
         ! in case dabs(pi*dt)>1.0d-6
       ELSE
         ! in case dabs(pi*dt)<=1.0d-6
         expp=expm1(-pi*dt)
         del2n=del2+del2*expp-expp*dyn2/pi
         del2int=-del2*expp/pi+dyn2*dt/pi+dyn2*expp/(pi*pi)
         del1n=+ri*del2*expp/pi-ri*dyn2*dt/pi- &
                    ri*dyn2*expp/(pi*pi)+dyn1*dt+del1
         del1int=-ri*del2*expp/(pi*pi)-ri*dyn2*dt*dt/(2.0d0*pi)+ &
                      ri*dyn2*expp/(pi*pi*pi)+dyn1*dt*dt/2.0d0+ &
                      del1*dt-ri*del2*dt/pi+ri*dyn2*dt/(pi*pi)
goto     100
       END IF
       ! ... only ice                                                      (end)
       ! in case rw==0.0d0
     END IF

     IF (irw == 1 .and. iri == 1) THEN
       a=(rw-pi)*(rw-pi)+4.0e0*pw*ri
       IF (a < 0.0) THEN
         PRINT*,   'in SUBROUTINE jersupsat: a < 0'
         WRITE (*,'(A,D13.5)')  'deter',        deter
         WRITE (*,'(A,4D13.5)') 'rw,pw,ri,pi',  rw,pw,ri,pi
         WRITE (*,'(A,3D13.5)') 'dt,dyn1,dyn2', dt,dyn1,dyn2
         WRITE (*,'(A,2D13.5)') 'del1,del2',    del1,del2
         PRINT*,   'stop 1905:a < 0'
         CALL finish(TRIM(modname),"fatal error: stop 1905:a < 0, model stop")
       END IF
       ! ... water and ice                                               (start)
       alfa=dsqrt((rw-pi)*(rw-pi)+4.0d0*pw*ri)
       ! beta is negative to the simple solution so it will decay

       beta=0.5d0*(alfa+rw+pi)
       gama=0.5d0*(alfa-rw-pi)
       g31=pi*dyn1-ri*dyn2
       g32=-pw*dyn1+rw*dyn2
       g2=rw*pi-ri*pw
       IF (g2 < 1.0d-20) g2 = 1.0004d-11*1.0003d-11-1.0002d-11*1.0001e-11 ! ... (ks) - 24th,may,2016
       expb=EXP(-beta*dt)
       expg=EXP(gama*dt)

       IF (dabs(gama*dt)>1.0e-6) THEN
         c11=(beta*del1-rw*del1-ri*del2-beta*g31/g2+dyn1)/alfa
         c21=(gama*del1+rw*del1+ri*del2-gama*g31/g2-dyn1)/alfa
         c12=(beta*del2-pw*del1-pi*del2-beta*g32/g2+dyn2)/alfa
         c22=(gama*del2+pw*del1+pi*del2-gama*g32/g2-dyn2)/alfa
         del1n=c11*expg+c21*expb+g31/g2
         del1int=c11*expg/gama-c21*expb/beta+(c21/beta-c11/gama) &
                 +g31*dt/g2
         del2n=c12*expg+c22*expb+g32/g2
         del2int=c12*expg/gama-c22*expb/beta+(c22/beta-c12/gama) &
                  +g32*dt/g2
goto     100
         ! in case dabs(gama*dt)>1.0d-6
       ELSE
         ! in case dabs(gama*dt)<=1.0d-6
         IF (abs(ri/rw)>1.0e-12) THEN
           IF (abs(rw/ri)>1.0e-12) THEN
             alfa=dsqrt((rw-pi)*(rw-pi)+4.0d0*pw*ri)
             beta=0.5d0*(alfa+rw+pi)
             gama=0.5d0*(alfa-rw-pi)
             IF (gama < 0.5*2.0d-10) gama=0.5d0*(2.002d-10-2.001d-10) ! ... (ks) - 24th,may,2016
             expg=expm1(gama*dt)
             expb=EXP(-beta*dt)
             ! beta/alfa could be very close to 1 that why i transform it
             ! remember alfa-beta=gama
             c11=(beta*del1-rw*del1-ri*del2+dyn1)/alfa
             c21=(gama*del1+rw*del1+ri*del2-gama*g31/g2-dyn1)/alfa
             c12=(beta*del2-pw*del1-pi*del2+dyn2)/alfa
             c22=(gama*del2+pw*del1+pi*del2-gama*g32/g2-dyn2)/alfa

             a1del1n=c11
             a2del1n=c11*expg
             a3del1n=c21*expb
             a4del1n=g31/g2*(gama/alfa+(gama/alfa-1.0d0)*expg)

             del1n=a1del1n+a2del1n+a3del1n+a4del1n

             a1del1int=c11*expg/gama
             a2del1int=-c21*expb/beta
             a3del1int=c21/beta
             a4del1int=g31/g2*dt*(gama/alfa)

             del1int=a1del1int+a2del1int+a3del1int+a4del1int

             a1del2n=c12
             a2del2n=c12*expg
             a3del2n=c22*expb
             a4del2n=g32/g2*(gama/alfa+ &
                        (gama/alfa-1.0d0)* &
                        (gama*dt+gama*gama*dt*dt/2.0d0))

             del2n=a1del2n+a2del2n+a3del2n+a4del2n

             a1del2int=c12*expg/gama
             a2del2int=-c22*expb/beta
             a3del2int=c22/beta
             a4del2int=g32/g2*dt*(gama/alfa)
             a5del2int=g32/g2*(gama/alfa-1.0d0)* &
                               (gama*dt*dt/2.0d0)

             del2int=a1del2int+a2del2int+a3del2int+a4del2int+ &
                     a5del2int
             ! in case dabs(rw/ri)>1d-12
           ELSE
             ! in case dabs(rw/ri)<=1d-12
             x=-2.0d0*rw*pi+rw*rw+4.0d0*pw*ri

             alfa=pi*(1+(x/pi)/2.0d0-(x/pi)*(x/pi)/8.0d0)
             beta=pi+(x/pi)/4.0d0-(x/pi)*(x/pi)/16.0d0+rw/2.0d0
             gama=(x/pi)/4.0d0-(x/pi)*(x/pi)/16.0d0-rw/2.0d0

             expg=expm1(gama*dt)
             expb=EXP(-beta*dt)

             c11=(beta*del1-rw*del1-ri*del2+dyn1)/alfa
             c21=(gama*del1+rw*del1+ri*del2-gama*g31/g2-dyn1)/alfa
             c12=(beta*del2-pw*del1-pi*del2+dyn2)/alfa
             c22=(gama*del2+pw*del1+pi*del2-gama*g32/g2-dyn2)/alfa

             del1n=c11+c11*expg+c21*expb+ &
                       g31/g2*(gama/alfa+(gama/alfa-1)*expg)
             del1int=c11*expg/gama-c21*expb/beta+(c21/beta)+ &
                         g31/g2*dt*(gama/alfa)
             del2n=c12+c12*expg+c22*expb+g32/g2*(gama/alfa+ &
                     (gama/alfa-1.0d0)* &
                     (gama*dt+gama*gama*dt*dt/2.0d0))
             del2int=c12*expg/gama-c22*expb/beta+ &
               (c22/beta)+g32/g2*dt*(gama/alfa)+ &
               g32/g2*(gama/alfa-1.0d0)*(gama*dt*dt/2.0d0)
             ! in case dabs(rw/ri)<=1d-12
           END IF
           ! alfa/beta 2
           ! in case dabs(ri/rw)>1d-12
         ELSE
           ! in case dabs(ri/rw)<=1d-12
           x=-2.0d0*rw*pi+pi*pi+4.0d0*pw*ri

           alfa=rw*(1.0d0+(x/rw)/2.0d0-(x/rw)*(x/rw)/8.0d0)
           beta=rw+(x/rw)/4.0d0-(x/rw)*(x/rw)/16.0d0+pi/2.0d0
           gama=(x/rw)/4.0d0-(x/rw)*(x/rw)/16.0d0-pi/2.0d0

           expg=expm1(gama*dt)
           expb=EXP(-beta*dt)

           c11=(beta*del1-rw*del1-ri*del2+dyn1)/alfa
           c21=(gama*del1+rw*del1+ri*del2-gama*g31/g2-dyn1)/alfa
           c12=(beta*del2-pw*del1-pi*del2+dyn2)/alfa
           c22=(gama*del2+pw*del1+pi*del2-gama*g32/g2-dyn2)/alfa

           del1n=c11+c11*expg+c21*expb+ &
                    g31/g2*(gama/alfa+(gama/alfa-1.0d0)*expg)
           del1int=c11*expg/gama-c21*expb/beta+(c21/beta)+ &
                      g31/g2*dt*(gama/alfa)
           del2n=c12+c12*expg+c22*expb+g32/g2* &
                    (gama/alfa+ &
                    (gama/alfa-1.0d0)*(gama*dt+gama*gama*dt*dt/2.0d0))
           del2int=c12*expg/gama-c22*expb/beta+c22/beta+ &
                g32/g2*dt*(gama/alfa)+ &
                g32/g2*(gama/alfa-1.0d0)*(gama*dt*dt/2.0d0)
           ! alfa/beta
           ! in case dabs(ri/rw)<=1d-12
         END IF
         ! in case dabs(gama*dt)<=1d-6
       END IF
       ! water and ice                                                 (end)
       ! in case isym1/=0.and.isym2/=0
     END IF

100  CONTINUE

     RETURN
   END SUBROUTINE jersupsat_ks
        
   SUBROUTINE jerdfun_ks (xi,xin,b21_my, &
                          fi2,psi2,del2n, &
                          isym2,ind,itype,idrop, &
                          nkr,col,ihydro,iin,kin)
     IMPLICIT NONE
     INTEGER,INTENT(IN) :: isym2, ind, itype, nkr, ihydro, iin, kin 
     INTEGER,INTENT(INOUT) :: idrop
     REAL(KIND=wp),INTENT(IN) :: b21_my(:), fi2(:), del2n, col
                                          ! fr_lim(:), frh_lim(:), &
     REAL(KIND=wp),INTENT(IN) :: xi(:)    ! tpn
     REAL(KIND=wp),INTENT(INOUT) :: xin(:)
     REAL(KIND=wp),INTENT(INOUT) :: psi2(:)
     INTEGER :: kr, nr, idsd_negative ! ityp
     REAL(KIND=wp) :: fi2r(nkr), psi2r(nkr), c, degree1, degree2, degree3, d, ratexi, &
                        b, a, xir(nkr),xinr(nkr) !, fr_lim_kr

     c = 2.0d0/3.0d0
     degree1 = 1.0d0/3.0d0
     degree2 = c
     degree3 = 3.0d0/2.0d0

!    IF (ind > 1) THEN
!      ityp = itype
!    ELSE
!      ityp = 1
!    END IF

     DO kr=1,nkr
       psi2r(kr) = fi2(kr)
       fi2r(kr) = fi2(kr)
     END DO
     nr=nkr
     ! new size distribution functions                             (start)
     IF (isym2 == 1) THEN
       IF (ind==1 .and. itype==1) THEN
         ! drop diffusional growth
         DO kr=1,nkr
           d=xi(kr)**degree1
           ratexi=c*del2n*b21_my(kr)/d
           b=xi(kr)**degree2
           a=b+ratexi
           IF (a<0.0d0) THEN
             xin(kr)=1.0d-50
           ELSE
             xin(kr)=a**degree3
           END IF
         END DO
         ! in case ind==1.and.itype==1
       ELSE
         ! in case ind/=1.OR.itype/=1
         DO kr=1,nkr
           ratexi = del2n*b21_my(kr)
           xin(kr) = xi(kr) + ratexi
         END DO
       END IF

       ! recalculation of size distribution FUNCTIONs                (start)
       DO kr=1,nkr
         xir(kr) = xi(kr)
         xinr(kr) = xin(kr)
         fi2r(kr) = fi2(kr)
       END DO

       idsd_negative = 0
       CALL jernewf_ks(nr,xir,fi2r,psi2r,xinr,isign_3point,idrop,nkr,col,idsd_negative,ihydro,iin,kin)
       IF (idsd_negative == 1)THEN
         IF (isign_ko_1 == 1) THEN
          ! we do not use kovatch-ouland as separate method
          ! CALL jernewf_ko_ks(nr,xir,fi2r,psi2r,xinr,nkr,col)
         END IF
       END IF

       DO kr=1,nkr
!        IF (itype==5) THEN
!          fr_lim_kr=frh_lim(kr)
!        ELSE
!          fr_lim_kr=fr_lim(kr)
!        END IF
         IF (psi2r(kr)<0.0d0) THEN
           PRINT*,    'stop 1506 : psi2r(kr)<0.0d0, in jerdfun_ks'
           CALL finish(TRIM(modname),"fatal error in psi2r(kr)<0.0d0, in jerdfun_ks, model stop")
         END IF
         psi2(kr) = psi2r(kr)
       END DO
       ! cycle by ice
       ! recalculation of size distribution FUNCTIONs                  (end)
       ! in case isym2/=0
     END IF
     ! new size distribution FUNCTIONs                               (end)

     RETURN
   END SUBROUTINE jerdfun_ks

   SUBROUTINE jernewf_ks(nrx,rr,fi,psi,rn,i3point,idrop,nkr,col,idsd_negative,ihydro, &
              iin,kin)
     IMPLICIT NONE
     INTEGER,INTENT(IN) :: nrx, i3point, nkr, ihydro, iin, kin
     INTEGER,INTENT(INOUT) :: idrop, idsd_negative
!del:    REAL(KIND=wp),INTENT(IN) :: tpn
     REAL(KIND=wp),INTENT(IN) :: col
     REAL(KIND=wp),INTENT(INOUT) :: psi(:), rn(:), fi(:), rr(:)
     INTEGER :: kmax, kr, i, k , nrxp, isign_diffusional_growth, nrx1,  &
              i3point_condevap !, ievap
     REAL(KIND=wp) :: rntmp,rrtmp,rrp,rrm,rntmp2,rrtmp2,rrp2,rrm2, gn1,gn2, &
                gn3,gn1p,gmat,gmat2, &
                cdrop(nrx),delta_cdrop(nrx),rrs(nrx+1),psinew(nrx+1), &
                psi_im,psi_i,psi_ip

     INTEGER,PARAMETER :: krdrop_remaping_min = 6, krdrop_remaping_max = 12

     nrxp = nrx + 1
     nrx1 = nkr

     DO i=1,nrx
       ! rn(i), g - new masses after condensation or evaporation
       IF (rn(i) < 0.0d0) THEN
         rn(i) = 1.0d-50
         fi(i) = 0.0d0
       END IF
     END DO

     DO k=1,nrx
       rrs(k)=rr(k)
     END DO

     i3point_condevap = i3point

!    ievap = 0
     IF (rn(1) < rrs(1)) THEN
       ! evaporation
       i3point_condevap = 0
       idrop = 0
       nrx1 = nrx
!      ievap = 1
     END IF

     IF (idrop == 0) i3point_condevap = 0

     DO k=1,nrx
       psi(k)=0.0d0
       cdrop(k)=0.0d0
       delta_cdrop(k)=0.0d0
       psinew(k)=0.0d0
     END DO
     rrs(nrxp)=rrs(nrx)*1024.0d0
     psinew(nrxp) = 0.0d0

     isign_diffusional_growth = 0
     DO k=1,nrx
       IF (rn(k).NE.rr(k)) THEN
       isign_diffusional_growth = 1
       goto 2000
       END IF
     END DO

2000 CONTINUE

     IF (isign_diffusional_growth == 1) THEN
     ! kovetz-olund method                                         (start)
       DO k=1,nrx1 ! nrx1-1
         IF (fi(k) > 0.0) THEN
           IF (dabs(rn(k)-rr(k)) < 1.0d-16) THEN
             psinew(k) = fi(k)*rr(k)
             cycle
           END IF

           i = 1
           DO WHILE (.not.(rrs(i) <= rn(k) .and. rrs(i+1) >= rn(k)) &
                  .and.i.LT.nrx1) ! was nrx1-1
             i = i + 1
           END DO

           IF (rn(k).LT.rrs(1)) THEN
             rntmp=rn(k)
             rrtmp=0.0d0
             rrp=rrs(1)
             gmat2=(rntmp-rrtmp)/(rrp-rrtmp)
             psinew(1)=psinew(1)+fi(k)*rr(k)*gmat2
           ELSE
             rntmp=rn(k)
             rrtmp=rrs(i)
             rrp=rrs(i+1)
             gmat2=(rntmp-rrtmp)/(rrp-rrtmp)
             gmat=(rrp-rntmp)/(rrp-rrtmp)
             psinew(i)=psinew(i)+fi(k)*rr(k)*gmat
             psinew(i+1)=psinew(i+1)+fi(k)*rr(k)*gmat2
           END IF
         END IF
!!3000     CONTINUE
       END DO

       DO kr=1,nrx1
         psi(kr)=psinew(kr)
       END DO

       DO kr=nrx1+1,nrx
         psi(kr)=fi(kr)
       END DO
       ! kovetz-olund method                                 (end)

       ! calculation both new total drop concentrations(after ko) and new
       ! total drop masses (after ko)

       ! 3point method	                                         (start)
       IF (i3point_condevap == 1) THEN
         DO k=1,nrx1-1
           IF (fi(k) > 0.0) THEN
             IF (dabs(rn(k)-rr(k)).LT.1.0d-16) THEN
               psi(k) = fi(k)*rr(k)
               goto 3001
             END IF

             IF (rrs(2).LT.rn(k)) THEN
               i = 2
               DO WHILE &
                 (.not.(rrs(i) <= rn(k) .and. rrs(i+1) >= rn(k)) &
                      .and.i.LT.nrx1-1)
                 i = i + 1
               END DO
               rntmp=rn(k)
               rrtmp=rrs(i)
               rrp=rrs(i+1)
               rrm=rrs(i-1)
               rntmp2=rn(k+1)
               rrtmp2=rrs(i+1)
               rrp2=rrs(i+2)
               rrm2=rrs(i)
               gn1=(rrp-rntmp)*(rrtmp-rntmp)/(rrp-rrm)/ &
                  (rrtmp-rrm)
               gn1p=(rrp2-rntmp2)*(rrtmp2-rntmp2)/ &
                    (rrp2-rrm2)/(rrtmp2-rrm2)
               gn2=(rrp-rntmp)*(rntmp-rrm)/(rrp-rrtmp)/ &
                    (rrtmp-rrm)
               gmat=(rrp-rntmp)/(rrp-rrtmp)
               gn3=(rrtmp-rntmp)*(rrm-rntmp)/(rrp-rrm)/ &
                                            (rrp-rrtmp)
               gmat2=(rntmp-rrtmp)/(rrp-rrtmp)
               psi_im = psi(i-1)+gn1*fi(k)*rr(k)
               psi_i = psi(i)+gn1p*fi(k+1)*rr(k+1)+&
                    (gn2-gmat)*fi(k)*rr(k)
               psi_ip = psi(i+1)+(gn3-gmat2)*fi(k)*rr(k)
               IF (psi_im > 0.0d0) THEN
                 IF (psi_ip > 0.0d0) THEN
                   IF (i > 2) THEN
                     ! smoothing criteria
                     IF (psi_im > psi(i-2) .and. psi_im < psi_i &
                       .and. psi(i-2) < psi(i) .OR. psi(i-2) >= psi(i)) THEN
                       psi(i-1) = psi_im
                       psi(i) = psi(i) + fi(k)*rr(k)*(gn2-gmat)
                       psi(i+1) = psi_ip
                     END IF
                   END IF
                 ELSE
                   exit
                 END IF
               ELSE
                 exit
               END IF
             END IF
           END IF
           
3001       CONTINUE

         END DO
       END IF
       ! 3 point method                                    (end)

       ! psi(k) - new hydrometeor size distribution FUNCTION
       DO k=1,nrx1
         psi(k)=psi(k)/rr(k)
       END DO
       DO k=nrx1+1,nrx
         psi(k)=fi(k)
       END DO
       IF (idrop == 1) THEN
         DO k=krdrop_remaping_min,krdrop_remaping_max
           cdrop(k)=3.0d0*col*psi(k)*rr(k)
         END DO
         ! kmax - right boundary spectrum of drop sdf
         !(krdrop_remap_min =< kmax =< krdrop_remap_max)
         DO k=krdrop_remaping_max,krdrop_remaping_min,-1
           kmax=k
           IF (psi(k).GT.0.0d0) goto 2011
         END DO

2011     CONTINUE
         DO k=kmax-1,krdrop_remaping_min,-1
           IF (cdrop(k).GT.0.0d0) THEN
             delta_cdrop(k)=cdrop(k+1)/cdrop(k)
             IF (delta_cdrop(k).LT.coeff_remaping) THEN
               cdrop(k)=cdrop(k)+cdrop(k+1)
               cdrop(k+1)=0.0d0
             END IF
           END IF
         END DO

         DO k=krdrop_remaping_min,kmax
           psi(k)=cdrop(k)/(3.0d0*col*rr(k))
         END DO
       END IF
       ! in case isign_diffusional_growth.NE.0
     ELSE
       ! in case isign_diffusional_growth.EQ.0
       DO k=1,nrx
         psi(k)=fi(k)
       END DO
     END IF

     DO kr=1,nrx
       IF (psi(kr) < 0.0) THEN ! ... (ks)
         idsd_negative = 1
         PRINT*, "idsd_negative=",idsd_negative,"kr",kr
         PRINT*,    'in SUBROUTINE jernewf'
         PRINT*,    'psi(kr)<0'
         PRINT*,    'before exit'
         PRINT*,    'isign_diffusional_growth'
         PRINT*,     isign_diffusional_growth
         PRINT*,    'i3point_condevap'
         PRINT*,     i3point_condevap
         PRINT*,    'k,rr(k),rn(k),k=1,nrx'
         PRINT*,    (k,rr(k),rn(k),k=1,nrx)
         PRINT*,    'k,rr(k),rn(k),fi(k),psi(k),k=1,nrx'
         WRITE (*,'(A,1X,I2,2X,4D13.5)') (k,rr(k),rn(k),fi(k),psi(k),k=1,nrx)
         PRINT*,    idrop,ihydro,iin,kin
         CALL finish(TRIM(modname),"fatal error in SUBROUTINE jernewf psi(kr)<0, < min, model stop")
       END IF
     END DO

     RETURN
   END SUBROUTINE jernewf_ks

   SUBROUTINE jerdfun_new_ks(xi,xin,b21_my,fi2,psi2, &
              idrop,nkr,col,ihydro,iin,kin)
     IMPLICIT NONE
     INTEGER,INTENT(INOUT) :: idrop, nkr
     INTEGER,INTENT(IN) :: ihydro,iin,kin
     REAL(KIND=wp),INTENT(IN) :: fi2(:), b21_my(:), col
     REAL(KIND=wp),INTENT(IN) :: xi(:)   ! tpn
     REAL(KIND=wp),INTENT(INOUT) :: psi2(:)
     REAL(KIND=wp),INTENT(INOUT) :: xin(:)
     INTEGER :: nr, kr, idsd_negative
     REAL(KIND=wp) :: c, degree1, degree2, degree3, d, ratexi, b, a, &
                        xir(nkr),fi2r(nkr),psi2r(nkr),xinr(nkr)
     c=2.0d0/3.0d0
     degree1=c/2.0d0
     degree2=c
     degree3=3.0d0/2.0d0
     nr=nkr
     xir = xi
     fi2r = fi2
     psi2r = psi2
     xinr = xin

     ! new drop size distribution functions                             (start)

     ! drop diffusional growth
     DO kr=1,nkr
       d = xir(kr)**degree1
       ratexi = c*b21_my(kr)/d
       b = xir(kr)**degree2
       a = b+ratexi
       IF (a<0.0d0) THEN
         xinr(kr) = 1.0d-50
       ELSE
         xinr(kr) = a**degree3
       END IF
     END DO

     ! recalculation of size distribution functions                (start)
     idsd_negative = 0
     !calculate the new PSD (remapping:
     CALL jernewf_ks(nr,xir,fi2r,psi2r,xinr,isign_3point,idrop,nkr,col,idsd_negative,ihydro,iin,kin)
     IF (idsd_negative == 1)THEN
       IF (isign_ko_2 == 1) THEN
       ! we do not use kovetz-ouland as separate method
       ! CALL jernewf_ko_ks(nr,xir,fi2r,psi2r,xinr,nkr,col)
       END IF
     END IF

     psi2 = psi2r
     ! recalculation of drop size distribution FUNCTIONs             (end)
     ! new drop size distribution FUNCTIONs                          (end)


     RETURN
   END SUBROUTINE jerdfun_new_ks

   SUBROUTINE jernucl01_ks(psi1_r,fccnr_r,fccnr_nucl_r,   &
                          xl_r,tt,ror_r,pp_r,sup1,&
                          col_r,rccn_r,nkr,nkr_aerosol,  &
                          win_r,is_this_cloudbase,ro_solute,ions,mwaero)
     IMPLICIT NONE
     INTEGER,INTENT(IN) :: nkr,nkr_aerosol,ions,is_this_cloudbase
     REAL(KIND=wp),INTENT(IN) ::  xl_r(:),ror_r,pp_r,col_r,rccn_r(:)
     REAL(KIND=wp),INTENT(IN) ::  mwaero,ro_solute,win_r
     REAL(KIND=wp),INTENT(INOUT) :: psi1_r(:),fccnr_r(:),fccnr_nucl_r(:)
     REAL(KIND=wp),INTENT(INOUT) :: tt,sup1
     REAL(KIND=wp) :: dropconcn(nkr), tpc, win
     REAL(KIND=wp) :: psi1(nkr),fccnr(nkr_aerosol),fccnr_nucl(nkr_aerosol),ror,xl(nkr),pp,col, &
             rccn(nkr_aerosol) !,dropradii(nkr)

!    oper3(ar1,ar2) = ar1*ar2/(0.622d0+0.378d0*ar1)


     ! ... adjust the input
     psi1 = psi1_r               !drop size distribution
     fccnr = fccnr_r             !ccn size distribution
     fccnr_nucl = fccnr_nucl_r   !nucleated ccn size distribution. PK: It is not used ?????
     xl = xl_r                   !water mass grid
     ror = ror_r                 !water density
     pp = pp_r                   !water pressure
     col = col_r                 !resolution of logaitmic mass grid: ln(2)/3
     rccn = rccn_r               !aerosol radii grid
!    dropradii = dropradii_r     !drop radii grid
     win = win_r

!    col3 = 3.0d0*col
!    rori = 1.0d0/ror

     ! -----------------------------
     ! ... drop nucleation (start)
     ! -----------------------------
!    tpn = tt
!    qpn = qq

     tpc = tt - 273.15d0

     IF (sup1>0.0d0 .and. tpc > t_nucl_drop_min) THEN !do nucleation if S over water >0 and t>t_nucl_drop_min=-80 celcius
       IF (sum(fccnr) > 0.0)THEN                      !if ccn concentration >0
         dropconcn = 0.0d0
         CALL water_nucleation (col, nkr_aerosol, psi1, fccnr, fccnr_nucl, xl, tt, ror, sup1, dropconcn, &
              pp, is_this_cloudbase, win, ro_solute, rccn, ions,mwaero)
       END IF
     END IF

  ! -----------------------------
  ! ... drop nucleation (end)
  ! -----------------------------


     ! ... output
     psi1_r = psi1              !new psd of water
     fccnr_r = fccnr            !new psd of ccn
     fccnr_nucl_r = fccnr_nucl  !new psd of nucleated ccn

     RETURN
   END SUBROUTINE jernucl01_ks

   SUBROUTINE water_nucleation (col, nkr, psi1, fccnr, fccnr_nucl, xl, tt, ror, sup1,     &
                            dropconcn, pp, is_this_cloudbase, win, ro_solute, &
                            rccn, ions, mwaero)
   !===================================================================!
   !                                                                   !
   ! drop nucleation scheme                                            !
   !                                                                   !
   ! authors: khain a.p. & pokrovsky a.g. july 2002 at huji, israel    !
   !                                                                   !
   !===================================================================!
     IMPLICIT NONE

     ! psi1(kr), 1/g/cm3 - non conservative drop size distribution function
     ! fccnr(kr), 1/cm^3 - aerosol(ccn) non conservative, size distribution function
     ! xl((kr), g        - drop bin masses

     INTEGER,INTENT(IN) :: is_this_cloudbase, nkr, ions
     REAL(KIND=wp),INTENT(IN) :: xl(:), ror, pp, win, rccn(:), col
     REAL(KIND=wp),INTENT(INOUT) :: fccnr(:), fccnr_nucl(:), psi1(:), dropconcn(:), tt, sup1
     REAL(KIND=wp),INTENT(IN) ::  ro_solute, mwaero
     INTEGER :: imax, i, ncriti, kr
     REAL(KIND=wp) :: dx,rcriti,deg01,rori,ccnconc(nkr),akoe,bkoe, rccn_minimum, &
               dln1, dln2, rmassl_nucl

!    oper3(ar1,ar2)=ar1*ar2/(0.622d0+0.378d0*ar1)
     dropconcn(:) = 0.0d0
     deg01 = 1.0d0/3.0d0
     rori=1.0/ror

     ! imax - right ccn spectrum boundary
     imax = nkr
     DO i=imax,1,-1
       IF (fccnr(i) > 0.0d0) THEN
         imax = i
         exit
       END IF
     END DO

     ncriti=0
     ! every iteration we will nucleate one bin, then we will check the new supersaturation
     ! and new rcriti.
     DO WHILE (imax>=ncriti)
       ccnconc = 0.0
       ! akoe & bkoe - constants in koehler equation
       akoe=3.3d-05/tt
       bkoe = ions*4.3/mwaero
       bkoe=bkoe*(4.0d0/3.0d0)*3.141593d0*ro_solute

       IF (use_cloud_base_nuc == 1) THEN
         IF (is_this_cloudbase == 1) THEN
           CALL cloud_base_super (fccnr, rccn, tt, pp, win, nkr, rcriti, ro_solute, ions, mwaero, col)
         ELSE
           ! rcriti, cm - critical radius of "dry" aerosol
           rcriti = (akoe/3.0d0)*(4.0d0/bkoe/sup1/sup1)**deg01
         END IF
       ELSE ! ismax_cloud_base==0
         ! rcriti, cm - critical radius of "dry" aerosol
         rcriti=(akoe/3.0d0)*(4.0d0/bkoe/sup1/sup1)**deg01
       END IF

       IF (rcriti >= rccn(imax)) exit ! nothing to nucleate
       ! find the minimum bin to nucleate
       ncriti = imax
       DO WHILE (rcriti<=rccn(ncriti) .and. ncriti>1)
         ncriti=ncriti-1
       END DO

       ! rccn_minimum - minimum aerosol(ccn) radius
       rccn_minimum = rccn(1)/10000.0d0
       ! calculation of ccnconc(ii)=fccnr(ii)*col - aerosol(ccn) bin concentrations, ii=imin,...,imax
       ! determination of ncriti   - number bin in which is located rcriti
       ! calculation of ccnconc(ncriti)=fccnr(ncriti)*dln1/(dln1+dln2),
       ! where,
       ! dln1=ln(rcriti)-ln(rccn_minimum)
       ! dln2=ln(rccn(1)-ln(rcriti)
       ! calculation of new value of fccnr(ncriti)

       ! the problem is that ncriti=1 and imax=2 may occur 2 times:
       ! a. rccn(1)<rcriti<rccn(2) --> we should clean part of the "1-2" bin: fccnr(imax)=fccnr(imax)*dln1/col
       ! b. rcriti<rccn(1)         --> we should clean the entire bin:  fccnr(imax)=0 
       ! and there is an option when criti=1 and imax=1, then:
       ! c. rcriti<rccn(1)         --> we should clean part of the "0-1" bin: fccnr(imax)=fccnr(imax)*dln1/(dln1+dln2)
       !                               this bin is special and has a width of dln1+dln2 and not just col
       !---------------------------------------------------
       IF ((imax-1>ncriti) .OR. ((ncriti==1) .and. (imax==2) .and. (rcriti<rccn(1)))) THEN ! imax bin should be cleaned to zero 
         ccnconc(imax) = col*fccnr(imax)
         fccnr_nucl(imax) = fccnr_nucl(imax) + fccnr(imax)
         fccnr(imax) = 0.0d0
       ELSE IF ((ncriti==1) .and. (imax==1) .and. (rcriti<rccn(1))) THEN ! rcriti<rccn(1) we should clean part of the "0-1" bin rccn_minimum<-->rccn(1)
         dln1=DLOG(rcriti/rccn_minimum)
         dln2=DLOG(rccn(1)/rcriti)
         ccnconc(imax)=dln2*fccnr(imax)
         fccnr_nucl(imax) = fccnr_nucl(imax) + fccnr(imax)*(1.0 - (dln1/DLOG(rccn(1)/rccn_minimum)))
         fccnr(imax)=fccnr(imax)*dln1/DLOG(rccn(1)/rccn_minimum)
       ELSE IF ((ncriti==imax-1) .and. (rcriti > rccn(imax-1))) THEN ! rccn(ncriti)<rcriti<rccn(imax) imax bin should be cleaned partially
         dln1=DLOG(rcriti/rccn(imax-1))
         dln2=col-dln1
         ccnconc(imax)=dln2*fccnr(imax)
         fccnr_nucl(imax) = fccnr_nucl(imax) + fccnr(imax)*(1.0 - dln1/col)
         fccnr(imax)=fccnr(imax)*dln1/col
       ELSE
         WRITE (txt,*) 'rcriti,rccn1,ncriti,imax=',rcriti,rccn(1),ncriti,imax
         CALL finish(TRIM(modname),'ccn bins problem, model stop' )
       END IF

       ! calculate the mass change due to nucleation
       rmassl_nucl=0.0d0
       IF (imax <= nkr-7) THEN ! we pass it to drops mass grid
         dropconcn(1) = dropconcn(1)+ccnconc(imax)
         rmassl_nucl = rmassl_nucl+ccnconc(imax)*xl(1)*xl(1)
       ELSE
         dropconcn(8-(nkr-imax)) = dropconcn(8-(nkr-imax))+ccnconc(imax)
         rmassl_nucl = rmassl_nucl + ccnconc(imax)*xl(8-(nkr-imax))*xl(8-(nkr-imax))
       END IF
       rmassl_nucl = rmassl_nucl*col*3.0*rori
       ! prepering to check if we need to nucleate the next bin
       imax = imax-1
       ! cycle imax>=ncriti
     END DO

     ! intergarting for including the previous nucleated drops
     IF (sum(dropconcn) > 0.0)THEN
       DO kr = 1,8
         dx = 3.0d0*col*xl(kr)
         psi1(kr) = psi1(kr)+dropconcn(kr)/dx
       END DO
     END IF

     RETURN
   END SUBROUTINE water_nucleation

   SUBROUTINE cloud_base_super (fccnr,rccn,tt,pp,wbase,nkr,rcriti,ro_solute,ions,mwaero,col)

     IMPLICIT NONE
     ! rccn(nkr),  cm- aerosol's radius
     ! fccnr(kr), 1/cm^3 - aerosol(ccn) non conservative, size
     !                     distribution function in point with x,z
     !                     coordinates, kr=1,...,nkr
     INTEGER,INTENT(IN) ::   nkr, ions
     REAL(KIND=wp),INTENT(IN) ::  tt, pp, wbase, rccn(:), col
     REAL(KIND=wp),INTENT(INOUT) :: fccnr(:), rcriti
     REAL(KIND=wp),INTENT(IN) ::  mwaero, ro_solute
     INTEGER :: nr, nn, kr
     REAL(KIND=wp) :: pl(nkr), supmax(nkr), akoe, bkoe, c3, pr, ccnconact, dl1, dl2

     CALL supmax_coeff(akoe,bkoe,c3,pp,tt,ro_solute,ions,mwaero)
     ! supmax calculation
     ! 'analytical estimation of droplet concentration at cloud base', eq.21, 2012
     ! calculation of right side hand of equation for s_max
     ! WHILE wbase>0, calculation pr

     pr = c3*wbase**(0.75d0)

     ! calculation supersaturation in cloud base
     supmax = 999.0
     pl = 0.0
     nn = -1
     DO nr=2,nkr
       supmax(nr)=dsqrt((4.0d0*akoe**3.0d0)/(27.0d0*rccn(nr)**3.0d0*bkoe))
       ! calculation ccnconact- the concentration of ccn that were activated
       ! following nucleation
       ! ccnconact=n from the paper
       ! 'analytical estimation of droplet concentration at cloud base', eq.19, 2012
       ! ccnconact, 1/cm^3- concentration of activated ccn = new droplet concentration
       ! ccnconact=fccnr(kr)*col
       ! col=ln2/3

       ccnconact=0.0d0

       ! nr represents the number of bin in which rcriti is located
       ! from nr bin to nkr bin goes to droplets
       DO kr=nr,nkr
         ccnconact = ccnconact + col*fccnr(kr)
       END DO
       ! calculate lhs of equation for s_max
       ! when pl<pr ccn will activate
       pl(nr)=supmax(nr)*(dsqrt(ccnconact))
       IF (pl(nr).LE.pr) THEN
         nn = nr
         exit
       END IF
     END DO ! nr

     IF (nn == -1) THEN
       PRINT*,"pr, wbase [cm/s], c3",pr,wbase,c3
       PRINT*,"pl",pl
       CALL finish(TRIM(modname),'nn is not defined in cloud base routine, model stop' )
     END IF

     ! linear interpolation- finding radius criti of aerosol between
     ! bin number (nn-1) to (nn)
     ! 1) finding the difference between pl and pr in the left and right over the
     ! final bin.

     dl1 = dabs(pl(nn-1)-pr) ! left side in the final bin
     dl2 = dabs(pl(nn)-pr)   ! right side in the final bin

     ! 2) fining the left part of bin that will not activate
     !	  dln1=col*dl1/(dl2+dl1)
     ! 3) finding the right part of bin that activate
     !	  dln2=col-dln1
     ! 4) finding radius criti of aerosol- rcriti

     rcriti = rccn(nn-1)*EXP(col*dl1/(dl1+dl2))
     ! end linear interpolation

     RETURN
   END SUBROUTINE cloud_base_super

   SUBROUTINE supmax_coeff (akoe,bkoe,c3,pp,tt,ro_solute,ions,mwaero)
     IMPLICIT NONE

     ! akoe, cm- constant in koehler equation
     ! bkoe    - constant in koehler equation
     ! f, cm^-2*s-  from koehler equation
     ! c3 - coefficient depends on thermodynamical parameters
     ! pp, (dynes/cm/cm)- pressure
     ! tt, (k)- temperature

     INTEGER,INTENT(IN) :: ions
     REAL(KIND=wp),INTENT(IN) :: pp, tt
     REAL(KIND=wp),INTENT(OUT) :: akoe, bkoe, c3
     REAL(KIND=wp),INTENT(IN) :: mwaero, ro_solute
     REAL(KIND=wp) , PARAMETER :: rv_my = 461.5d4, cp = 1005.0d4, g = 9.8d2, rd_my = 287.0d4, & ![cgs]
                          pi = 3.141593d0
     REAL(KIND=wp) :: pzero,tzero,alw1,sw,ro_w,hc,ew,ro_v,dv,ro_a,fl,fr,f,tpc,qv,a1,a2, &
               c1 !,c2,deg01,deg02

     tpc = tt-273.15d0
     ! cgs :
     ! rv_my, cm*cm/sec/sec/kelvin - individual gas constant
     !                               for water vapour
     ! rv_my=461.5d4
     ! cp,  cm*cm/sec/sec/kelvin- specific heat capacity of
     !	                            moist air at constant pressure
     ! cp=1005.0d4
     ! g,  cm/sec/sec- acceleration of gravity
     ! g=9.8d2
     ! rd_my, cm*cm/sec/sec/kelvin - individual gas constant
     !                               for dry air
     ! rd_my=287.0d4
     ! al2_my, cm*cm/sec/sec - latent heat of sublimation
     ! al2_my=2.834d10
     ! pzero, dynes/cm/cm - reference pressure
     pzero=1.01325d6
     ! tzero, kelvin - reference temperature
     tzero=273.15d0
     ! al1_my, cm*cm/sec/sec - latent heat of vaporization
     ! alw1=al1_my - alw1 depends on temperature
     ! alw1, [m^2/sec^2] -latent heat of vaporization-
     alw1 = -6.143419998d-2*tpc**(3.0d0)+1.58927d0*tpc**(2.0d0) &
          -2.36418d3*tpc+2.50079d6
     ! alw1, [cm^2/sec^2]
     alw1 = alw1*10.0d3
     ! sw, [n*m^-1] - surface tension of water-air interface
     IF (tpc.LT.-5.5d0) THEN
       sw = 5.285d-11*tpc**(6.0d0)+6.283d-9*tpc**(5.0d0)+ &
          2.933d-7*tpc**(4.0d0)+6.511d-6*tpc**(3.0d0)+ &
          6.818d-5*tpc**(2.0d0)+1.15d-4*tpc+7.593d-2
     ELSE
       sw = -1.55d-4*tpc+7.566165d-2
     END IF
     ! sw, [g/sec^2]
     sw = sw*10.0d2
     ! ro_w, [kg/m^3] - density of liquid water
     IF (tpc.LT.0.0d0) THEN
       ro_w= -7.497d-9*tpc**(6.0d0)-3.6449d-7*tpc**(5.0d0) &
           -6.9987d-6*tpc**(4.0d0)+1.518d-4*tpc**(3.0d0) &
           -8.486d-3*tpc**(2.0d0)+6.69d-2*tpc+9.9986d2
     ELSE
       ro_w=(-3.932952d-10*tpc**(5.0d0)+1.497562d-7*tpc**(4.0d0) &
           -5.544846d-5*tpc**(3.0d0)-7.92221d-3*tpc**(2.0d0)+ &
           1.8224944d1*tpc+9.998396d2)/(1.0d0+1.8159725d-2*tpc)
     END IF
     ! ro_w, [g/cm^3]
     ro_w=ro_w*1.0d-3
     ! hc, [kg*m/kelvin*sec^3] - coefficient of air heat conductivity
     hc=7.1128d-5*tpc+2.380696d-2
     ! hc, [g*cm/kelvin*sec^3]
     hc=hc*10.0d4

     ! ew, water vapor pressure ! ... ks (kg/m2/sec)
     ew = 6.38780966d-9*tpc**(6.0d0)+2.03886313d-6*tpc**(5.0d0)+ &
        3.02246994d-4*tpc**(4.0d0)+2.65027242d-2*tpc**(3.0d0)+ &
        1.43053301d0*tpc**(2.0d0)+4.43986062d1*tpc+6.1117675d2

     ! ew, [g/cm*sec^2]
     ew=ew*10.0d0
     ! akoe & bkoe - constants in koehler equation
     !ro_solute=2.16d0
     akoe=2.0d0*sw/(rv_my*ro_w*(tpc+tzero))
     bkoe = ions*4.3/mwaero
     bkoe=bkoe*(4.0d0/3.0d0)*pi*ro_solute

     ! ro_v, g/cm^3 - density of water vapor
     !                calculate from equation of state for water vapor
     ro_v = ew/(rv_my*(tpc+tzero))
     ! dv,  [cm^2/sec] - coefficient of diffusion

     ! 'pruppacher, h.r., klett, j.d., 1997. microphysics of clouds and precipitation'
     ! 'page num 503, eq. 13-3'
     dv = 0.211d0*(pzero/pp)*((tpc+tzero)/tzero)**(1.94d0)

     ! qv,  g/g- water vapor mixing ratio
     ! ro_a, g/cm^3 - density of air, from equation of state
     ro_a=pzero/((tpc+tzero)*rd_my)

     ! f, s/m^2 - coefficient depending on thermodynamics parameters
     !            such as temperature, thermal conductivity of air, etc
     ! left side of f equation
     fl=(ro_w*alw1**(2.0d0))/(hc*rv_my*(tpc+tzero)**(2.0d0))

     ! right side of f equation
     fr = ro_w*rv_my*(tpc+tzero)/(ew*dv)
     f = fl + fr

     ! qv, g/g - water vapor mixing ratio
     qv=ro_v/ro_a

     ! a1,a2 -  terms from equation describing changes of
     !          supersaturation in an adiabatic cloud air parcel
     ! a1, [cm^-1] - constant
     ! a2, [-]     - constant

     a1=(g*alw1/(cp*rv_my*(tpc+tzero)**(2.0d0)))-(g/(rd_my*(tpc+tzero)))
     a2=(1.0d0/qv)+(alw1**(2.0d0))/(cp*rv_my*(tpc+tzero)**(2.0d0))

     ! c1,c2,c3,c4- constant parameters

     c1=1.058d0
!    c2=1.904d0
!    deg01=1.0d0/3.0d0
!    deg02=1.0d0/6.0d0
     c3=c1*(f*a1/3.0d0)**(0.75d0)*dsqrt(3.0d0*ro_a/(4.0d0*pi*ro_w*a2))

     RETURN
   END SUBROUTINE supmax_coeff

   SUBROUTINE coll_xxx_bott(g,ckxx,x,c,ima,prdkrn,nkr,output_flux)

     IMPLICIT NONE
     INTEGER,INTENT(IN) :: nkr
     REAL(KIND=wp),INTENT(INOUT) :: g(:)
     REAL(KIND=wp),INTENT(IN) :: ckxx(:,:),x(:), c(:,:)
     INTEGER,INTENT(IN) :: ima(:,:)
     REAL(KIND=wp),INTENT(IN) :: prdkrn
     REAL(KIND=wp),INTENT(INOUT) :: output_flux(:)
     REAL(KIND=wp):: gmin,x01,x02,x03,gsi,gsj,gsk,gk,flux,x1
     INTEGER :: i,ix0,ix1,j,k,kp
!    INTEGER :: kp_flux_max

     gmin=1.0d-16
!    kp_flux_max = nkr

     ! ix0 - lower limit of integration by i
     DO i=1,nkr-1
       ix0=i
       IF (g(i).GT.gmin) goto 2000
     END DO
2000 CONTINUE
     IF (ix0.EQ.nkr-1) RETURN

     ! ix1 - upper limit of integration by i
     DO i=nkr-1,1,-1
       ix1=i
       IF (g(i).GT.gmin) goto 2010
     END DO
2010 CONTINUE

     IF (ix1 == nkr) ix1 = nkr - 1

     ! ... collisions
     DO i=ix0,ix1
       IF (g(i).LE.gmin) goto 2020
       DO j=i,ix1
         IF (g(j).LE.gmin) goto 2021
         k=ima(i,j)
         kp=k+1
         x01=ckxx(i,j)*g(i)*g(j)*prdkrn
         x02=dmin1(x01,g(i)*x(j))
         IF (j.NE.k) x03=dmin1(x02,g(j)*x(i))
         IF (j.EQ.k) x03=x02
         gsi=x03/x(j)
         gsj=x03/x(i)
         gsk=gsi+gsj
         g(i)=g(i)-gsi  ! this needs to be limited (for all the hydro)
         g(j)=g(j)-gsj
         gk=g(k)+gsk    ! when j=/k needs to be limited (only for different hydro)
         flux=0.d0

         ! g(i), g(j) - PSD f of bins i,j.
         ! their parts sum up (=gsk) have to be added to bin gk.
         ! since it falls between k and kp=k+1, it is added to bins g(k) and g(kp=k+1) via type of remapping (see flux below)

         IF (gk.GT.gmin) THEN
           x1=DLOG(g(kp)/gk+1.0d-12) ! avoid LOG(1) --> x1=0
           flux=gsk/x1*(EXP(0.5*x1)-EXP(x1*(0.5-c(i,j))))
           flux=min(flux,gk,gsk)
           g(k)=gk-flux
           g(kp)=g(kp)+flux
           ! --- [js] output flux - for autoconv.
           output_flux(kp) = output_flux(kp) + flux
         END IF

         IF (g(i) < 0.0 .OR. g(j) < 0.0 .OR. g(k) < 0.0 .OR. g(kp) < 0.0) THEN 
           PRINT*,    'i,j,k,kp'
           PRINT*,     i,j,k,kp
           PRINT*,    'ix0,ix1'
           PRINT*,     ix0,ix1
           PRINT*,   'g(i),g(j),g(k),g(kp)'
           WRITE(*,'(A,4D13.5)') g(i),g(j),g(k),g(kp)
           stop 'stop in collisions'
         END IF
2021     CONTINUE
       END DO
2020   CONTINUE
     END DO

     RETURN
   END SUBROUTINE coll_xxx_bott

   SUBROUTINE coll_xxx_bott_mod1(g,is,ie,je,ckxx,x,c,ima,prdkrn,nkr,output_flux)
     IMPLICIT NONE
     INTEGER,INTENT(IN) :: nkr,is,ie,je
     REAL(KIND=wp),INTENT(INOUT) :: g(:)
     REAL(KIND=wp),INTENT(IN) :: ckxx(:,:),x(:), c(:,:)
     INTEGER,INTENT(IN) :: ima(:,:)
     REAL(KIND=wp),INTENT(IN) :: prdkrn
     REAL(KIND=wp),INTENT(INOUT) :: output_flux(:)
     REAL(KIND=wp):: gmin,x01,x02,x03,gsi,gsj,gsk,gk, flux,x1
     INTEGER :: i,ix0,ix1,j,k,kp,iee,jee
!    INTEGER :: kp_flux_max

     gmin=1.0d-16
!    kp_flux_max = nkr

     ! ix0 - lower limit of integration by i
     DO i=1,nkr-1
       ix0=i
       IF (g(i).GT.gmin) goto 2000
     END DO
2000 CONTINUE
     IF (ix0.EQ.nkr-1) RETURN
     
     ! ix1 - upper limit of integration by i
     DO i=nkr-1,1,-1
       ix1=i
       IF (g(i).GT.gmin) goto 2010
     END DO
2010 CONTINUE

     iee = ie; jee = je
     IF (iee == nkr) iee = nkr-1
     IF (jee == nkr) jee = nkr-1

     ! ... collisions
     DO i=is,iee
       IF (g(i).LE.gmin) goto 2020
       DO j=i,jee
         IF (g(j).LE.gmin) goto 2021
         k=ima(i,j)
         kp=k+1
         x01=ckxx(i,j)*g(i)*g(j)*prdkrn
         x02=dmin1(x01,g(i)*x(j))
         IF (j.NE.k) x03=dmin1(x02,g(j)*x(i))
         IF (j.EQ.k) x03=x02
         gsi=x03/x(j)
         gsj=x03/x(i)
         gsk=gsi+gsj
         g(i)=g(i)-gsi
         g(j)=g(j)-gsj
         gk=g(k)+gsk
         flux=0.d0
         ! g(i), g(j) - psd f of bins i,j. 
         ! their parts sum up (=gsk) have to be added to bin gk. 
         ! since it falls between k and kp=k+1, it is added to bins g(k) and g(kp=k+1) via type of remapping (see flux below)
         IF (gk.GT.gmin) THEN
           x1=DLOG(g(kp)/gk+1.0d-12) ! avoid LOG(1) --> x1=0
           flux=gsk/x1*(EXP(0.5*x1)-EXP(x1*(0.5-c(i,j))))
           flux=min(flux,gk,gsk)
           g(k)=gk-flux
           g(kp)=g(kp)+flux
           ! --- [js] output flux - for autoconv.
           output_flux(kp) = output_flux(kp) + flux
         END IF

         IF (g(i) < 0.0 .OR. g(j) < 0.0 .OR. g(k) < 0.0 .OR. g(kp) < 0.0) THEN 
           PRINT*,    'i,j,k,kp'
           PRINT*,     i,j,k,kp
           PRINT*,    'ix0,ix1'
           PRINT*,     ix0,ix1
           PRINT*,   'g(i),g(j),g(k),g(kp)'
           WRITE(*,'(A,4D13.5)') g(i),g(j),g(k),g(kp)
           stop 'stop in collisions'
         END IF
2021     CONTINUE
       END DO
2020   CONTINUE
     END DO

     RETURN
   END SUBROUTINE coll_xxx_bott_mod1



   SUBROUTINE coll_xxx_bott_mod2(g,is,ie,js,je,ckxx,x,c,ima,prdkrn,nkr,output_flux)
     IMPLICIT NONE
     INTEGER,INTENT(IN) :: nkr,is,js,ie,je
     REAL(KIND=wp),INTENT(INOUT) :: g(:)
     REAL(KIND=wp),INTENT(IN) :: ckxx(:,:),x(:), c(:,:)
     INTEGER,INTENT(IN) :: ima(:,:)
     REAL(KIND=wp),INTENT(IN) :: prdkrn
     REAL(KIND=wp),INTENT(INOUT) :: output_flux(:)
     REAL(KIND=wp):: gmin,x01,x02,x03,gsi,gsj,gsk,gk, flux,x1
     INTEGER :: i,ix0,ix1,j,k,kp,iee,jee
!    INTEGER :: kp_flux_max
  
     gmin=1.0d-16
!    kp_flux_max = nkr
  
     ! ix0 - lower limit of integration by i
     DO i=1,nkr-1
       ix0=i
       IF (g(i).GT.gmin) goto 2000
     END DO
2000 CONTINUE
     IF (ix0.EQ.nkr-1) RETURN
  
     ! ix1 - upper limit of integration by i
     DO i=nkr-1,1,-1
       ix1=i
       IF (g(i).GT.gmin) goto 2010
     END DO
2010 CONTINUE
  
     iee = ie; jee = je
     IF (iee == nkr) iee = nkr-1
     IF (jee == nkr) jee = nkr-1

     ! ... collisions
     DO i=is,iee
       IF (g(i).LE.gmin) goto 2020
       DO j=js,jee
         IF (g(j).LE.gmin) goto 2021
         k=ima(i,j)
         kp=k+1
         x01=ckxx(i,j)*g(i)*g(j)*prdkrn
         x02=dmin1(x01,g(i)*x(j))
         IF (j.NE.k) x03=dmin1(x02,g(j)*x(i))
         IF (j.EQ.k) x03=x02
         gsi=x03/x(j)
         gsj=x03/x(i)
         gsk=gsi+gsj
         g(i)=g(i)-gsi
         g(j)=g(j)-gsj
         gk=g(k)+gsk
         flux=0.d0
         ! g(i), g(j) - psd f of bins i,j.
         ! their parts sum up (=gsk) have to be added to bin gk.
         ! since it falls between k and kp=k+1, it is added to bins g(k) and g(kp=k+1) via type of remapping (see flux below)

         IF (gk.GT.gmin) THEN
           x1=DLOG(g(kp)/gk+1.0d-12) ! avoid LOG(1) --> x1=0
           flux=gsk/x1*(EXP(0.5*x1)-EXP(x1*(0.5-c(i,j))))
           flux=min(flux,gk,gsk)
           g(k)=gk-flux
           g(kp)=g(kp)+flux
           ! --- [js] output flux - for autoconv.
           output_flux(kp) = output_flux(kp) + flux
         END IF
  
         IF (g(i) < 0.0 .OR. g(j) < 0.0 .OR. g(k) < 0.0 .OR. g(kp) < 0.0) THEN  
           PRINT*,    'i,j,k,kp'
           PRINT*,     i,j,k,kp
           PRINT*,    'ix0,ix1'
           PRINT*,     ix0,ix1
           PRINT*,   'g(i),g(j),g(k),g(kp)'
           WRITE(*,'(A,4D13.5)') g(i),g(j),g(k),g(kp)
           stop 'stop in collisions'
         END IF
2021     CONTINUE
       END DO
2020   CONTINUE
     END DO

     RETURN
   END SUBROUTINE coll_xxx_bott_mod2
END MODULE mo_sbm_main 
