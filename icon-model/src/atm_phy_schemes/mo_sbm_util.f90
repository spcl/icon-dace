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

MODULE mo_sbm_util

  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: finish, message, txt => message_text
  USE mo_mpi,                ONLY: my_process_is_stdio, p_bcast, p_comm_work, p_io
  USE mo_io_units,           ONLY: find_next_free_unit
  USE mo_nwp_tuning_config,  ONLY: tune_sbmccn

  USE mo_model_domain,       ONLY: t_patch
  USE mo_nonhydro_types,     ONLY: t_nh_prog
  USE mo_run_config,         ONLY: iqb_s, msg_level
  USE mo_impl_constants,     ONLY: min_rlcell
  USE mo_loopindices,        ONLY: get_indices_c
  USE mo_atm_phy_nwp_config, ONLY: atm_phy_nwp_config
  USE mo_2mom_mcrph_driver,  ONLY: two_moment_mcrph_init

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: &
       & warm_hucminit,              &
       & kernals_ks,                 &
       & ccn_init_sbm,               &
       & sbm_init,                   &
       & ibreakup, snow_breakup_on, spont_rain_breakup_on, &
       & conserv, jiwen_fan_melt, ipolar_hucm, hail_opt,   &
       & dx_bound, scal, iceprocs, iceturb, liqturb,       &
       & icempl, icemax, ncd, nhydr, nhydro, k0_ll, krmin_ll, &
       & l0_ll, ieps_400, ieps_800, ieps_1600, k0l_gl, k0g_gl,&
       & krminl_gl, krmaxl_gl, krming_gl, krmaxg_gl, kr_icempl,&
       & krbreak, krice, krdrop, nkr, jmax, nrg, jbreak, br_max, &
       & krmin_breakup, nkr_aerosol, dt_coll, c1_mey, c2_mey, &
       & col, p1, p2, p3, alcr, alcr_g, ncond, ncoll, kp_flux_max, &
       & g_lim, kr_sgs_max, isign_ko_1, isign_ko_2, isign_3point, &
       & idebug_print_debugmodule, &
       & coeff_remaping, ventpl_max, rw_pw_min, ri_pi_min, rw_pw_ri_pi_min, ratio_icew_min, &
       & use_cloud_base_nuc, t_nucl_drop_min, t_nucl_ice_min, ice_nucl_method, isign_tq_icenucl, &
       & delsupice_max, &
       & radxx, massxx, denxx, massxxo, denxxo, vri, xx, roccn, fccnr_mix, fccnr, ff1r_d, xl_d, vr1_d, &
       & ff3r_d,xs_d,vr3_d,vts_d,fliqfr_sd,ro3bl_d, &
       & ff4r_d,xg_d,vr4_d,vtg_d,fliqfr_gd,ro4bl_d, &
       & ff5r_d,xh_d,vr5_d,vth_d,fliqfr_hd,ro5bl_d, &
       & xs_melt_d,xg_melt_d,xh_melt_d,vr_test,frimfr_sd,rf3r, &
       & xi_melt_d, ff2r_d, xi_d, vr2_d, vtc_d, fliqfr_id, ro2bl_d, &
       & t_new_d, rhocgs_d, pcgs_d, dt_d, qv_old_d, qv_d, &
       & xl_mg, &
       & bin_mass, tab_colum, tab_dendr, tab_snow, bin_log, &
       & rlec, rsec, rgec, rhec, xl, xs, xg, xh, vr1, vr3, vr4, vr5, &
       & riec, xi, vr2, coefin, slic, tlic, ywll_1000mb, ywll_750mb, ywll_500mb, &
       & ro1bl, ro2bl, ro3bl, ro4bl, ro5bl, radxxo, &
       & ima, chucm, brkweight, ecoalmassm, prob, gain_var_new, nnd, &
       & dropradii, pkij, qkj, ikr_spon_break, cwll, &
       & fccnr_mar, fccnr_con, fccnr_obs, ccnr, scale_ccn_factor, xccn, &
       & rccn, &
       & icloud, ilognormal_modes_aerosol, do_aero_bc, iccn_reg, mwaero, &
       & ions, ro_solute, diagccn, fccnorig, fccnd, &
       & fr_lim, frh_lim,lh_ce_1, lh_ce_2, lh_ce_3,  &
       & lh_frz, lh_mlt, lh_rime, lh_homo, ce_bf, ce_af, ds_bf, &
       & ds_af, mlt_bf, mlt_af, frz_af, frz_bf, cldnucl_af,     &
       & cldnucl_bf, icenucl_af, icenucl_bf, lh_ice_nucl, del_cldnucl_sum, &
       & del_icenucl_sum, del_ce_sum, del_ds_sum,del_ccnreg, ccn_reg, &
       & auto_cld_nsink_b,auto_cld_msink_b, &
       & accr_cld_nsink_b,accr_cld_msink_b,  &
       & selfc_rain_nchng_b,dbl_orhocgs,dbl_odt, &
       & ttdiffl, automass_ch, autonum_ch, nrautonum, &
       & p_ff1i01, p_ff1i33, p_ff8i01, p_ff8i33, p_ff8in01, p_ff8in33, &
       & z0in, zmin

 CHARACTER(LEN=*), PARAMETER :: modname = 'mo_sbm_util'
 INTEGER, PARAMETER ::          p_ff1i01=1,   p_ff1i33=33, & !droplets
                                p_ff8i01=34,  p_ff8i33=66, & !ccn
                                p_ff8in01=67, p_ff8in33=99   !nucleated ccn for regeneration
 INTEGER,PARAMETER :: ibreakup = 1
 INTEGER,PARAMETER :: snow_breakup_on = 1
 INTEGER,PARAMETER :: spont_rain_breakup_on = 1
 LOGICAL, PARAMETER :: conserv = .true.
 INTEGER,PARAMETER :: jiwen_fan_melt = 1
 LOGICAL, PARAMETER :: ipolar_hucm = .false.
 INTEGER,PARAMETER :: hail_opt = 1
 REAL, PARAMETER :: dx_bound = 1433
 REAL(KIND=wp), PARAMETER ::  scal = 1.d0
 INTEGER,PARAMETER :: iceprocs = 0
 INTEGER,PARAMETER :: iceturb = 0, liqturb = 0
 INTEGER,PARAMETER :: icempl=1,icemax=3,ncd=33,nhydr=5,nhydro=7                &
                ,k0_ll=8,krmin_ll=1,krmax_ll=19,l0_ll=6                  &
                ,ieps_400=1,ieps_800=0,ieps_1600=0                       &
                ,k0l_gl=16,k0g_gl=16                                     &
                ,krminl_gl=1,krmaxl_gl=24                                &
                ,krming_gl=1,krmaxg_gl=33,kr_icempl=9                    &
                ,krbreak=17,krice=18,krdrop=15                           & ! krdrop=bin 15 --> 50um
                !,nkr=43,jmax=43,nrg=2,jbreak=28,br_max=43,krmin_breakup=31,nkr_aerosol=43   ! option for 43 bins
                ,nkr=33,jmax=33,nrg=2,jbreak=18,br_max=33,krmin_breakup=31,nkr_aerosol=33    ! option for 33 bins
 REAL(KIND=wp), PARAMETER :: z0in=2.0e5,zmin=2.0e5
 REAL(KIND=wp) :: dt_coll
 REAL(KIND=wp), PARAMETER :: c1_mey=0.00033,c2_mey=0.0,col=0.23105, &
                   p1=1000000.0,p2=750000.0,p3=500000.0,  &
                   alcr = 0.5, &
                   alcr_g = 100.0 ! forcing no transition from graupel to hail in this version
 INTEGER :: ncond, ncoll

 INTEGER,PARAMETER :: kp_flux_max = 33
 REAL(KIND=wp), PARAMETER :: g_lim = 1.0d-16 ! [g/cm^3]
 INTEGER,PARAMETER :: kr_sgs_max = 20 ! rg(20)=218.88 mkm

 INTEGER,PARAMETER :: isign_ko_1 = 0, isign_ko_2 = 0,  isign_3point = 1,  &
                      idebug_print_debugmodule = 1
 DOUBLE PRECISION, PARAMETER::coeff_remaping = 0.0066667d0
 DOUBLE PRECISION, PARAMETER::ventpl_max = 5.0d0

 DOUBLE PRECISION, PARAMETER::rw_pw_min = 1.0d-10
 DOUBLE PRECISION, PARAMETER::ri_pi_min = 1.0d-10
 DOUBLE PRECISION, PARAMETER::rw_pw_ri_pi_min = 1.0d-10
 DOUBLE PRECISION, PARAMETER::ratio_icew_min = 1.0d-4

 INTEGER,PARAMETER :: use_cloud_base_nuc = 0
 REAL(KIND=wp), PARAMETER::t_nucl_drop_min = -80.0d0
 REAL(KIND=wp), PARAMETER::t_nucl_ice_min = -37.0d0
!ice nucleation method
!using meyers method  : ice_nucl_method == 0
!using de_mott method : ice_nucl_method == 1
 INTEGER,PARAMETER :: ice_nucl_method = 0
 INTEGER,PARAMETER :: isign_tq_icenucl = 1
 DOUBLE PRECISION, PARAMETER::delsupice_max = 59.0d0 !delsupice_max=59%
 REAL(KIND=wp) :: radxx(nkr,nhydr-1),massxx(nkr,nhydr-1),denxx(nkr,nhydr-1) &
                  ,massxxo(nkr,nhydro),denxxo(nkr,nhydro),vri(nkr)           &
                  ,xx(nkr),roccn(nkr),fccnr_mix(nkr),fccnr(nkr)
 REAL(KIND=wp),DIMENSION (nkr) :: ff1r_d,xl_d,vr1_d &
                            ,ff3r_d,xs_d,vr3_d,vts_d,fliqfr_sd,ro3bl_d &
                            ,ff4r_d,xg_d,vr4_d,vtg_d,fliqfr_gd,ro4bl_d &
                            ,ff5r_d,xh_d,vr5_d,vth_d,fliqfr_hd,ro5bl_d &
                            ,xs_melt_d,xg_melt_d,xh_melt_d,vr_test,frimfr_sd,rf3r
 ! ... radar variables
 REAL(KIND=wp),DIMENSION (nkr,icemax) :: xi_melt_d &
                            ,ff2r_d,xi_d,vr2_d,vtc_d,fliqfr_id,ro2bl_d
 REAL(KIND=wp) :: t_new_d,rhocgs_d,pcgs_d,dt_d,qv_old_d,qv_d
 REAL(KIND=wp) :: xl_mg(nkr)
 REAL(KIND=wp),PRIVATE :: & 
             & xs_mg(nkr),xg_mg(nkr),xh_mg(nkr), &
             & xi1_mg(nkr),xi2_mg(nkr),xi3_mg(nkr)
 ! ----------------------------------------------------------------------------------+
 ! ... warm-sbm-init
 ! ... holding lookup tables and memory arrays for the fast_sbm module
 REAL(KIND=wp), ALLOCATABLE, DIMENSION(:)::                             &
                                  bin_mass,tab_colum,tab_dendr,tab_snow,bin_log
 REAL(KIND=wp), ALLOCATABLE, DIMENSION(:) ::                            &
                                  rlec,rsec,rgec,rhec,xl,xs,xg,xh,vr1,vr3,vr4,vr5
 REAL(KIND=wp), ALLOCATABLE, DIMENSION(:,:)::                           &
                                  riec,xi,vr2
 REAL(KIND=wp), ALLOCATABLE ::                              &
                                  coefin(:),slic(:,:),tlic(:,:), &
                                  ywll_1000mb(:,:),ywll_750mb(:,:),ywll_500mb(:,:)
 REAL(KIND=wp), ALLOCATABLE ::                                                   &
                                  ro1bl(:), ro2bl(:,:), ro3bl(:), ro4bl(:), ro5bl(:),  &
                                  radxxo(:,:)
 INTEGER,ALLOCATABLE ::           ima(:,:)
 REAL(KIND=wp), ALLOCATABLE ::   chucm(:,:)
 REAL(KIND=wp), ALLOCATABLE ::   brkweight(:),ecoalmassm(:,:), prob(:),gain_var_new(:,:),nnd(:,:)
 REAL(KIND=wp), ALLOCATABLE ::   dropradii(:),pkij(:,:,:),qkj(:,:)
 INTEGER ::                       ikr_spon_break
 REAL(KIND=wp), ALLOCATABLE ::   cwll(:,:)
 REAL(KIND=wp),ALLOCATABLE ::     fccnr_mar(:),fccnr_con(:)
 REAL(KIND=wp),ALLOCATABLE ::     fccnr_obs(:),ccnr(:)
 REAL(KIND=wp),ALLOCATABLE ::     scale_ccn_factor,xccn(:),rccn(:)
 ! ... warm-sbm-init
 ! --------------------------------------------------------------------------------+  
 INTEGER :: icloud
! ----aerosol setup
! aerosol size distribution (sd)
 INTEGER,PARAMETER :: ilognormal_modes_aerosol = 1 !follow lognormal distribution
 ! in case of ilognormal_modes_aerosol = 0 ! read in a sd file from observation. Currently the file name for the observed
 ! sd is "ccn_size_33bin.dat", which is from the july 18 2017 "ena" case.
 INTEGER,PARAMETER :: do_aero_bc = 0
 INTEGER,PARAMETER :: iccn_reg = 1
 ! aerosol composition
 REAL(KIND=wp), PARAMETER :: mwaero = 22.9 + 35.5 ! sea salt. ! mwaero = 115.0
 INTEGER,PARAMETER :: ions = 2          ! sea salt. Set ions = 3 for ammonium-sulfate
 REAL(KIND=wp), PARAMETER :: ro_solute = 2.16        ! sea salt. Set ro_solute = 1.79 for ammonium-sulfate
 !for diagnostic ccn for places where sources exist
 LOGICAL, PARAMETER :: diagccn = .false.
 REAL(KIND=wp), ALLOCATABLE :: fccnorig(:), fccnd(:) ! for diagccn
 ! ----aerosol setup (end)

 REAL(KIND=wp) :: fr_lim(nkr), frh_lim(nkr),lh_ce_1, lh_ce_2, lh_ce_3,  &
                      lh_frz, lh_mlt, lh_rime, lh_homo, ce_bf, ce_af, ds_bf, &
                      ds_af, mlt_bf, mlt_af, frz_af, frz_bf, cldnucl_af,     &
                      cldnucl_bf, icenucl_af, icenucl_bf, lh_ice_nucl, del_cldnucl_sum, &
                      del_icenucl_sum, del_ce_sum, del_ds_sum,del_ccnreg, ccn_reg

 REAL(KIND=wp) :: auto_cld_nsink_b,auto_cld_msink_b, &
                       accr_cld_nsink_b,accr_cld_msink_b,  &
                       selfc_rain_nchng_b,dbl_orhocgs,dbl_odt

 !-----------------------------------------------------------
 REAL(KIND=wp) :: ttdiffl, automass_ch, autonum_ch, nrautonum

 ! following are two type declarations to hold the values
 ! of equidistand lookup tables. up to now, such lookup tables are
 ! used in the segal-khain parameterization of ccn-activation and
 ! for determining the wet growth diameter of graupel.

 CONTAINS

  SUBROUTINE warm_hucminit(fccnr_con,fccnr_mar,fccnr_obs)
    IMPLICIT NONE
    INTEGER :: unitnr
    INTEGER :: i,j,kr
    REAL(KIND=wp) :: dlnr,ax,deg01,concccnin,contccnin
    CHARACTER(LEN=256), PARAMETER :: dir_43 = "SBM_input_43", dir_33 = "SBM_input_33"
    CHARACTER(LEN=256) :: input_dir
    INTEGER :: error
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM(modname)//'::warm_hucminit'
    REAL(KIND=wp) ,INTENT(INOUT) :: fccnr_con(:), fccnr_mar(:), fccnr_obs(:)

    IF (nkr == 33) input_dir = TRIM(dir_33)
    IF (nkr == 43) input_dir = TRIM(dir_43)
    CALL message(routine," fast sbm: initializing hujisbm ")
    CALL message(routine," fast sbm: ****** hujisbm ******* ")

 ! +-------------------------------------------------------+
 ! lookuptable #1
 ! +-------------------------------------------------------+
    IF (.not. allocated(bin_mass)) allocate(bin_mass(nkr))
    IF (.not. allocated(tab_colum)) allocate(tab_colum(nkr))
    IF (.not. allocated(tab_dendr)) allocate(tab_dendr(nkr))
    IF (.not. allocated(tab_snow)) allocate(tab_snow(nkr))
    IF (.not. allocated(bin_log)) allocate(bin_log(nkr))

    dlnr=LOG(2.d0)/(3.d0)

    IF (my_process_is_stdio()) THEN
      unitnr=find_next_free_unit(10,999)
      WRITE(txt, '(a,i2)') 'warm_hucminit : table-1 -- opening BLKD_SDC.dat'
      CALL message(modname,TRIM(txt))
      OPEN(unitnr, file=TRIM(input_dir)//"/BLKD_SDC.dat", status='old', form='formatted', iostat=error)
      IF (error /= 0) THEN
        WRITE (txt,*) 'warm_hucminit : table-1 not found'
        CALL message(modname,TRIM(txt))
        CALL finish(TRIM(modname),'error in warm_hucminit')
      END IF
      DO kr=1,nkr
         READ (unitnr,*) bin_mass(kr),tab_colum(kr),tab_dendr(kr),tab_snow(kr)
         bin_log(kr) = LOG10(bin_mass(kr))
      END DO
    END IF

    CALL p_bcast(bin_mass, p_io, p_comm_work)
    CALL p_bcast(tab_colum, p_io, p_comm_work)
    CALL p_bcast(tab_dendr, p_io, p_comm_work)
    CALL p_bcast(tab_snow, p_io, p_comm_work)
    CALL p_bcast(bin_log, p_io, p_comm_work)

    WRITE(txt, '(a,i2)') 'fast_sbm_init : succesfull reading table-1'
    CALL message(routine,TRIM(txt))
 ! +-----------------------------------------------------------------------+
 ! lookuptable #2
 ! +----------------------------------------------+
    IF (.not. allocated(rlec)) allocate(rlec(nkr))
    IF (.not. allocated(riec)) allocate(riec(nkr,icemax))
    IF (.not. allocated(rsec)) allocate(rsec(nkr))
    IF (.not. allocated(rgec)) allocate(rgec(nkr))
    IF (.not. allocated(rhec)) allocate(rhec(nkr))

    IF (my_process_is_stdio()) THEN
      unitnr=find_next_free_unit(10,999)
      WRITE(txt, '(a,i2)') 'warm_hucminit : table-2 -- opening capacity33.asc'
      CALL message(modname,TRIM(txt))
      OPEN(unitnr, file=TRIM(input_dir)//"/capacity33.asc", status='old', form='formatted', iostat=error)
      IF (error /= 0) THEN
        WRITE (txt,*) 'warm_hucminit : table-2 not found'
        CALL message(modname,TRIM(txt))
        CALL finish(TRIM(modname),'error in warm_hucminit')
      END IF
      READ (unitnr,*) rlec,riec,rsec,rgec,rhec
    END IF

    CALL p_bcast(rlec, p_io, p_comm_work)
    CALL p_bcast(riec, p_io, p_comm_work)
    CALL p_bcast(rsec, p_io, p_comm_work)
    CALL p_bcast(rgec, p_io, p_comm_work)
    CALL p_bcast(rhec, p_io, p_comm_work)

    WRITE(txt, '(a,i2)') 'fast_sbm_init : succesfull reading table-2'
    CALL message(routine,TRIM(txt))
 ! +----------------------------------------------------------------------+
 ! lookuptable #3
 ! +-----------------------------------------------+
    IF (.not. allocated(xl)) allocate(xl(nkr))
    IF (.not. allocated(xi)) allocate(xi(nkr,icemax))
    IF (.not. allocated(xs)) allocate(xs(nkr))
    IF (.not. allocated(xg)) allocate(xg(nkr))
    IF (.not. allocated(xh)) allocate(xh(nkr))

    IF (my_process_is_stdio()) THEN
      unitnr=find_next_free_unit(10,999)
      WRITE(txt, '(a,i2)') 'warm_hucminit : table-3 -- opening masses33.asc'
      CALL message(modname,TRIM(txt))
      OPEN(unitnr, file=TRIM(input_dir)//"/masses33.asc", status='old', form='formatted', iostat=error)
      IF (error /= 0) THEN
        WRITE (txt,*) 'warm_hucminit : table-3 not found'
        CALL message(modname,TRIM(txt))
        CALL finish(TRIM(modname),'error in warm_hucminit')
      END IF
      READ (unitnr,*) xl,xi,xs,xg,xh
    END IF

    CALL p_bcast(xl, p_io, p_comm_work)
    CALL p_bcast(xi, p_io, p_comm_work)
    CALL p_bcast(xs, p_io, p_comm_work)
    CALL p_bcast(xg, p_io, p_comm_work)
    CALL p_bcast(xh, p_io, p_comm_work)

    WRITE(txt, '(a,i2)') 'fast_sbm_init : succesfull reading table-3'
    CALL message(routine,TRIM(txt))
 ! +-------------------------------------------------------------------------+
 ! lookuptable #4
 ! terminal velosity :
 ! +---------------------------------------------------+
    IF (.not. allocated(vr1)) allocate(vr1(nkr))
    IF (.not. allocated(vr2)) allocate(vr2(nkr,icemax))
    IF (.not. allocated(vr3)) allocate(vr3(nkr))
    IF (.not. allocated(vr4)) allocate(vr4(nkr))
    IF (.not. allocated(vr5)) allocate(vr5(nkr))

    IF (my_process_is_stdio()) THEN
      unitnr=find_next_free_unit(10,999)
      WRITE(txt, '(a,i2)') 'warm_hucminit : table-4 -- opening termvels.asc'
      CALL message(modname,TRIM(txt))
      OPEN(unitnr, file=TRIM(input_dir)//"/termvels33_corrected.asc", status='old', form='formatted', iostat=error)
      IF (error /= 0) THEN
        WRITE (txt,*) 'warm_hucminit : table-4 not found'
        CALL message(modname,TRIM(txt))
        CALL finish(TRIM(modname),'error in warm_hucminit')
      END IF
      READ (unitnr,*) vr1,vr2,vr3,vr4,vr5
    END IF

    CALL p_bcast(vr1, p_io, p_comm_work)
    CALL p_bcast(vr2, p_io, p_comm_work)
    CALL p_bcast(vr3, p_io, p_comm_work)
    CALL p_bcast(vr4, p_io, p_comm_work)
    CALL p_bcast(vr5, p_io, p_comm_work)

    WRITE(txt, '(a,i2)') 'fast_sbm_init : succesfull reading table-4'
    CALL message(routine,TRIM(txt))
 ! +----------------------------------------------------------------------+
 ! lookuptable #5
 ! constants : *** not need it here ***!
 ! +----------------------------------------------------------------------+

 ! +----------------------------------------------------------------------+
 ! lookuptable #6
 ! kernels depending on pressure :
 ! +------------------------------------------------------------------+
    IF (.not. allocated(ywll_1000mb)) allocate(ywll_1000mb(nkr,nkr))
    IF (.not. allocated(ywll_750mb)) allocate(ywll_750mb(nkr,nkr))
    IF (.not. allocated(ywll_500mb)) allocate(ywll_500mb(nkr,nkr))

    IF (my_process_is_stdio()) THEN
      unitnr=find_next_free_unit(10,999)
      WRITE(txt, '(a,i2)') 'warm_hucminit : table-6 -- opening kernels_z.asc'
      CALL message(modname,TRIM(txt))
      OPEN(unitnr, file=TRIM(input_dir)//"/kernLL_z33.asc", status='old', form='formatted', iostat=error)
      IF (error /= 0) THEN
        WRITE (txt,*) 'warm_hucminit : table-6 not found'
        CALL message(modname,TRIM(txt))
        CALL finish(TRIM(modname),'error in warm_hucminit')
      END IF
      READ (unitnr,*) ywll_1000mb,ywll_750mb,ywll_500mb
    END IF
    DO i=1,nkr
       DO j=1,nkr
          IF (i > 33 .or. j > 33) THEN
             ywll_1000mb(i,j) = 0.0
             ywll_750mb(i,j) =  0.0
             ywll_500mb(i,j) =  0.0
          END IF
       END DO
    END DO

    CALL p_bcast(ywll_1000mb, p_io, p_comm_work)
    CALL p_bcast(ywll_750mb, p_io, p_comm_work)
    CALL p_bcast(ywll_500mb, p_io, p_comm_work)

    WRITE(txt, '(a,i2)') 'fast_sbm_init : succesfull reading table-6'
    CALL message(routine,TRIM(txt))

 ! +-----------------------------------------------------------------------+
 ! lookuptable #7
 ! collisions kernels :  no ice-ice / ice-liquid here
 ! +-----------------------------------------------------------------------+

 ! +-----------------------------------------------------------------------+
 ! lookuptable #8
 ! bulkdensity:
 ! +--------------------------------------------------------------+
    IF (.not. allocated(ro1bl)) allocate(ro1bl(nkr))
    IF (.not. allocated(ro2bl)) allocate(ro2bl(nkr,icemax))
    IF (.not. allocated(ro3bl)) allocate(ro3bl(nkr))
    IF (.not. allocated(ro4bl)) allocate(ro4bl(nkr))
    IF (.not. allocated(ro5bl)) allocate(ro5bl(nkr))

    IF (my_process_is_stdio()) THEN
      unitnr=find_next_free_unit(10,999)
      WRITE(txt, '(a,i2)') 'warm_hucminit : table-8 -- opening bulkdens.asc'
      CALL message(modname,TRIM(txt))
      OPEN(unitnr, file=TRIM(input_dir)//"/bulkdens33.asc", status='old', form='formatted', iostat=error)
      IF (error /= 0) THEN
        WRITE (txt,*) 'warm_hucminit : table-8 not found'
        CALL message(modname,TRIM(txt))
        CALL finish(TRIM(modname),'error in warm_hucminit')
      END IF
      READ (unitnr,*) ro1bl,ro2bl,ro3bl,ro4bl,ro5bl
    END IF

    CALL p_bcast(ro1bl, p_io, p_comm_work)
    CALL p_bcast(ro2bl, p_io, p_comm_work)
    CALL p_bcast(ro3bl, p_io, p_comm_work)
    CALL p_bcast(ro4bl, p_io, p_comm_work)
    CALL p_bcast(ro5bl, p_io, p_comm_work)

    WRITE(txt, '(a,i2)') 'fast_sbm_init : succesfull reading table-8'
    CALL message(routine,TRIM(txt))
 ! +----------------------------------------------------------------------+
 ! lookuptable #9
 ! bulkradii:
 ! +-----------------------------------------------------------+
    IF (.not. allocated(radxxo)) allocate(radxxo(nkr,nhydro))

    IF (my_process_is_stdio()) THEN
      unitnr=find_next_free_unit(10,999)
      WRITE(txt, '(a,i2)') 'warm_hucminit : table-9 -- opening bulkradii.asc'
      CALL message(modname,TRIM(txt))
      OPEN(unitnr, file=TRIM(input_dir)//"/bulkradii33.asc", status='old', form='formatted', iostat=error)
      IF (error /= 0) THEN
        WRITE (txt,*) 'warm_hucminit : table-9 not found'
        CALL message(modname,TRIM(txt))
        CALL finish(TRIM(modname),'error in warm_hucminit')
      END IF
      READ (unitnr,*) radxxo
    END IF

    CALL p_bcast(radxxo, p_io, p_comm_work)

    WRITE(txt, '(a,i2)') 'fast_sbm_init : succesfull reading table-9'
    CALL message(routine,TRIM(txt))
 ! +-----------------------------------------------------------------------+
 ! lookuptable #10
 ! polar-hucm scattering amplitudes look-up table :
 ! +-----------------------------------------------------------------------+

 ! +-----------------------------------------------------------------------+

 ! calculation of the mass(in mg) for categories boundaries :
    ax=2.d0**(1.0)

    DO i=1,nkr
       xl_mg(i) = xl(i)*1.e3
       xs_mg(i) = xs(i)*1.e3
       xg_mg(i) = xg(i)*1.e3
       xh_mg(i) = xh(i)*1.e3
       xi1_mg(i) = xi(i,1)*1.e3
       xi2_mg(i) = xi(i,2)*1.e3
       xi3_mg(i) = xi(i,3)*1.e3
    END DO

    IF (.not. allocated(ima)) allocate(ima(nkr,nkr))
    IF (.not. allocated(chucm)) allocate(chucm(nkr,nkr))
    chucm  = 0.0d0
    ima = 0
    CALL courant_bott_ks(xl, nkr, chucm, ima, scal)
    WRITE(txt, '(a,i2)') 'fast_sbm_init : succesfull reading "courant_bott_ks" '
    CALL message(routine,TRIM(txt))

    deg01=1./3.
    concccnin=0.
    contccnin=0.
    IF (.not. allocated(dropradii)) allocate(dropradii(nkr))
    DO kr=1,nkr
       dropradii(kr)=(3.0*xl(kr)/4.0/3.141593/1.0)**deg01
    END DO

 ! +-------------------------------------------------------------+
 ! allocating aerosols array
 ! +-------------------------+
    IF (.not. allocated(xccn)) allocate(xccn(nkr_aerosol))
    IF (.not. allocated(rccn)) allocate(rccn(nkr_aerosol))
    IF (.not. allocated(scale_ccn_factor)) allocate(scale_ccn_factor)
    IF (.not. allocated(ccnr)) allocate(ccnr(nkr_aerosol))

! ... initializing the fccnr_mar and fccnr_con
    fccnr_con = 0.0
    fccnr_mar = 0.0
    fccnr_obs = 0.0
    scale_ccn_factor = 1.0
    xccn = 0.0
    rccn = 0.0

    IF (ilognormal_modes_aerosol == 1)THEN
       CALL lognormal_modes_aerosol(fccnr_con,fccnr_mar,nkr_aerosol,col,xl,xccn,rccn,scale_ccn_factor,1)
       CALL lognormal_modes_aerosol(fccnr_con,fccnr_mar,nkr_aerosol,col,xl,xccn,rccn,scale_ccn_factor,2)
       WRITE(txt, '(a,i2)') 'fast_sbm_init : succesfull reading "lognormal_modes_aerosol" '
       CALL message(routine,TRIM(txt))
!---yz2020mar:read aerosol size distribution from observation----
    ELSE ! read an observed sd with a format of aerosol size (cm), dn (#cm-3) and dndlogd for 33bins (jinwe fan)
       IF (my_process_is_stdio()) THEN
         unitnr=find_next_free_unit(10,999)
         WRITE(txt, '(a,i2)') 'warm_hucminit : opening CCN_size_33bin.dat'
         CALL message(modname,TRIM(txt))
         OPEN(unitnr, file=TRIM(input_dir)//"/CCN_size_33bin.dat", status='old', form='formatted', iostat=error)
         IF (error /= 0) THEN
           WRITE (txt,*) 'warm_hucminit : CCN_size_33bin.dat not found'
           CALL message(modname,TRIM(txt))
           CALL finish(TRIM(modname),'error in warm_hucminit')
         END IF
         DO kr=1,nkr
           READ (unitnr,*) rccn(kr),ccnr(kr),fccnr_obs(kr) !---aerosol size (cm), dn (# cm-3) and dndlogd for 33bins
         END DO
         CALL message(routine,"fast_sbm_init: succesfull reading aerosol sd from observation")
       END IF
    END IF

 ! +-------------------------------------------------------------+

    IF (.not. allocated(pkij)) allocate(pkij(jbreak,jbreak,jbreak))
    IF (.not. allocated(qkj)) allocate(qkj(jbreak,jbreak))
    IF (.not. allocated(ecoalmassm)) allocate(ecoalmassm(nkr,nkr))
    IF (.not. allocated(brkweight)) allocate(brkweight(jbreak))
    pkij = 0.0e0
    qkj = 0.0e0
    ecoalmassm = 0.0d0
    brkweight = 0.0d0
    CALL breakinit_ks(pkij,qkj,ecoalmassm,xl,dropradii,jbreak,nkr,vr1) ! rain spontanous breakup

    CALL p_bcast(pkij, p_io, p_comm_work)
    CALL p_bcast(qkj, p_io, p_comm_work)
    WRITE(txt, '(a,i2)') 'fast_sbm_init : succesfull reading breakinit_ks" '
    CALL message(routine,TRIM(txt))
  ! +--------------------------------------------------------------------------------------------------------------------+
    IF (.not. allocated(cwll)) allocate(cwll(nkr,nkr))
    cwll(:,:) = 0.0e0

!   100   FORMAT(10i4)
!   101   FORMAT(3x,f7.5,e13.5)
!   102   FORMAT(4e12.4)
!   105   FORMAT(a48)
!   106   FORMAT(a80)
!   123   FORMAT(3e12.4,3i4)
!   200   FORMAT(6e13.5)
!   201   FORMAT(6d13.5)
!   300   FORMAT(8e14.6)
!   301   FORMAT(3x,f8.3,3x,e13.5)
!   302   FORMAT(5e13.5)

    RETURN

    WRITE( txt , '(a,i4)' )                                          &
          'module_mp_fast_sbm_init: error opening hujisbm_data on unit,model stop '
    CALL finish(TRIM(modname),txt)

  END SUBROUTINE warm_hucminit

  SUBROUTINE kernals_ks(dtime_coal,nkr,p_z)
    IMPLICIT NONE
    INTEGER :: nkr
    REAL(KIND=wp),INTENT(IN) :: dtime_coal,p_z
    INTEGER :: i,j
    REAL(KIND=wp), PARAMETER :: p1=1.0e6,p2=0.75e6,p3=0.50e6
    REAL(KIND=wp) :: dlnr, scal, p_1, p_2, p_3, ckern_1, ckern_2, &
                                                  ckern_3
    scal = 1.0
    dlnr = LOG(2.0d0)/(3.0d0*scal)

    p_1=p1
    p_2=p2
    p_3=p3
    DO i=1,nkr
      DO j=1,nkr
        ! 1. water - water
        ckern_1 = ywll_1000mb(i,j)
        ckern_2 = ywll_750mb(i,j)
        ckern_3 = ywll_500mb(i,j)
        cwll(i,j) = ckern_z(p_z,p_1,p_2,p_3,ckern_1,ckern_2,ckern_3)*dtime_coal*dlnr
      END DO
    END DO

    ! ... ecoalmassm is from "breakiniit_ks"
    DO i=1,nkr
      DO j=1,nkr
        cwll(i,j) = ecoalmassm(i,j)*cwll(i,j)
      END DO
    END DO

    RETURN
  END SUBROUTINE kernals_ks

  REAL(KIND=wp) FUNCTION ckern_z (p_z,p_1,p_2,p_3,ckern_1,ckern_2,ckern_3)

    IMPLICIT NONE

    REAL(KIND=wp),INTENT(IN) :: p_z,p_1,p_2,p_3,ckern_1,ckern_2,ckern_3
    IF (p_z>=p_1) ckern_z = ckern_1
    !IF (p_z==p_2) ckern_z=ckern_2
    IF (p_z<=p_3) ckern_z = ckern_3
    IF (p_z<p_1 .and. p_z>=p_2) ckern_z = ckern_2 + (ckern_1-ckern_2)*(p_z-p_2)/(p_1-p_2)
    IF (p_z<p_2 .and. p_z>p_3) ckern_z = ckern_3 + (ckern_2-ckern_3)*(p_z-p_3)/(p_2-p_3)

    RETURN
  END FUNCTION ckern_z

  SUBROUTINE lognormal_modes_aerosol(fccnr_con,fccnr_mar,nkr_local,col,xl,xccn,rccn,scale_fa,itype)

    IMPLICIT NONE
    ! ... interface
    INTEGER,INTENT(IN) :: nkr_local, itype
    REAL(KIND=wp) ,INTENT(IN) :: xl(:), col, scale_fa
    REAL(KIND=wp) ,INTENT(OUT) :: fccnr_con(:), fccnr_mar(:)
    REAL(KIND=wp) ,INTENT(OUT) :: xccn(:),rccn(:)
    ! ... interface
    ! ... local
    INTEGER :: kr
    REAL(KIND=wp)  :: ccncon1, ccncon2, ccncon3, radius_mean1, radius_mean2, radius_mean3,  &
                      sig1, sig2, sig3
    REAL(KIND=wp)  :: concccnin, fccnr_tmp(nkr_local), deg01, x0drop,                           &
                      x0, r0, x0ccn, roccn(nkr_local),  &
                      arg11,arg12,arg13,arg21,arg22,arg23,      &
                      arg31,arg32,arg33,dnbydlogr_norm1,dnbydlogr_norm2,dnbydlogr_norm3

    REAL(KIND=wp) , PARAMETER :: rccn_max = 0.4d-4         ! [cm]
    ! ... minimal radii for dry aerosol for the 3 log normal distribution
    REAL(KIND=wp) , PARAMETER :: rccn_min_3ln = 0.00048d-4 ! [cm]
    REAL(KIND=wp) , PARAMETER :: pi = 3.14159265d0
    REAL(KIND=wp) , PARAMETER :: roccn0 = 0.1000e01 !---yz2020mar

    ! ... calculating the ccn radius grid
    !ro_solute_nacl = 2.16d0  ! [g/cm3]
    !ro_solute_ammon = 1.79  ! [g/cm3]

    ! note: rccn(1) = 1.2 nm
    !       rccn(33) = 2.1 um
    deg01 = 1.0d0/3.0d0
    x0drop = xl(1)
    x0ccn = x0drop/(2.0**(nkr_local))
    DO kr = nkr_local,1,-1
      roccn(kr) = roccn0
      x0 = x0ccn*2.0d0**(kr)
      r0 = (3.0d0*x0/4.0d0/3.141593d0/roccn(kr))**deg01
      xccn(kr) = x0
      rccn(kr) = r0
    END DO

    IF (itype == 1) THEN ! maritime regime

      ccncon1 = 340.000
      radius_mean1 = 0.00500d-04
      sig1 = 1.60000

      ccncon2 = 60.0000
      radius_mean2 = 0.03500d-04
      sig2 = 2.00000

      ccncon3 = 3.10000
      radius_mean3 = 0.31000d-04
      sig3 = 2.70000

    ELSE IF(itype == 2) THEN ! continental regime
      ccncon1 = 1000.000
!     IF (tune_sbmccn .LT. 0.06 ) THEN ! assuming clean case
!       ccncon1 = 6000.000 !then, when multiplied by tune_sbmccn=0.059, it will become ~340
!     END IF
      radius_mean1 = 0.00800d-04
      sig1 = 1.60000

      ccncon2 = 800.0000
      radius_mean2 = 0.03400d-04
      sig2 = 2.10000

      ccncon3 = 0.72000
      radius_mean3 = 0.46000d-04
      sig3 = 2.20000
    END IF

    fccnr_tmp = 0.0
    concccnin = 0.0

    arg11 = ccncon1/(sqrt(2.0d0*pi)*LOG(sig1))
    arg21 = ccncon2/(sqrt(2.0d0*pi)*LOG(sig2))
    arg31 = ccncon3/(sqrt(2.0d0*pi)*LOG(sig3))

    dnbydlogr_norm1 = 0.0
    dnbydlogr_norm2 = 0.0
    dnbydlogr_norm3 = 0.0
    DO kr = nkr_local,1,-1
      IF (rccn(kr) > rccn_min_3ln .and. rccn(kr) < rccn_max)THEN
        arg12 = (LOG(rccn(kr)/radius_mean1))**2.0
        arg13 = 2.0d0*((LOG(sig1))**2.0);
        dnbydlogr_norm1 = arg11*exp(-arg12/arg13)*(LOG(2.0)/3.0)
        arg22 = (LOG(rccn(kr)/radius_mean2))**2.0
        arg23 = 2.0d0*((LOG(sig2))**2.0)
        dnbydlogr_norm2 = dnbydlogr_norm1 + arg21*exp(-arg22/arg23)*(LOG(2.0)/3.0)
        arg32 = (LOG(rccn(kr)/radius_mean3))**2.0
        arg33 = 2.0d0*((LOG(sig3))**2.0)
        dnbydlogr_norm3 = dnbydlogr_norm2 + arg31*exp(-arg32/arg33)*(LOG(2.0)/3.0);
        fccnr_tmp(kr) = dnbydlogr_norm3
      END IF
    END DO
    concccnin = sum(fccnr_tmp(:))
    PRINT*,'concccnin',concccnin
    IF (itype == 1) fccnr_mar = scale_fa*fccnr_tmp/col
    IF (itype == 2) fccnr_con = scale_fa*fccnr_tmp/col

    RETURN
  END SUBROUTINE lognormal_modes_aerosol

     ! +----------------------------------------------------+
  SUBROUTINE courant_bott_ks(xl, nkr, chucm, ima, scal)

    IMPLICIT NONE

    INTEGER,INTENT(IN) :: nkr
    REAL(KIND=wp),INTENT(IN) :: xl(:)
    REAL(KIND=wp),INTENT(INOUT) :: chucm(:,:)
    INTEGER,INTENT(INOUT) :: ima(:,:)
    REAL(KIND=wp),INTENT(IN) :: scal

    ! ... locals
    INTEGER :: k, kk, j, i
    REAL(KIND=wp) :: x0, xl_mg(nkr), dlnr
    ! ... locals

    ! ima(i,j) - k-category number,
    ! chucm(i,j)   - courant number :
    ! logarithmic grid distance(dlnr) :

    !xl_mg(0)=xl_mg(1)/2
    xl_mg(1:nkr) = xl(1:nkr)*1.0d3

    dlnr=LOG(2.0d0)/(3.0d0*scal)

    DO i = 1,nkr
      DO j = i,nkr
        x0 = xl_mg(i) + xl_mg(j)
        DO k = j,nkr
          !IF (k == 1) goto 1000 ! ### (ks)
          kk = k
          IF (k == 1) goto 1000
            IF (xl_mg(k) >= x0 .and. xl_mg(k-1) < x0) THEN
              chucm(i,j) = LOG(x0/xl_mg(k-1))/(3.d0*dlnr)
              IF (chucm(i,j) > 1.0d0-1.d-08) THEN
                chucm(i,j) = 0.0d0
                kk = kk + 1
              END IF
              ima(i,j) = min(nkr-1,kk-1)
              !IF (ima(i,j) == 0) THEN
              !       PRINT*,"ima==0"
              !END IF
              goto 2000
            END IF
          1000 continue
        END DO
        2000  continue
        !IF (i.EQ.nkr.or.j.EQ.nkr) ima(i,j)=nkr
        chucm(j,i) = chucm(i,j)
        ima(j,i) = ima(i,j)
      END DO
    END DO

    RETURN
  END SUBROUTINE courant_bott_ks

  SUBROUTINE breakinit_ks(pkij,qkj,ecoalmassm,xl_r,dropradii,jbreak,nkr,vr1)
    IMPLICIT NONE
    ! ... interface
    INTEGER,INTENT(IN) :: jbreak, nkr
    REAL(KIND=wp),INTENT(INOUT) :: ecoalmassm(:,:)
    REAL(KIND=wp),INTENT(IN) :: xl_r(:), dropradii(:), vr1(:)
    REAL(KIND=wp),INTENT(INOUT) :: pkij(:,:,:),qkj(:,:)
    INTEGER :: unitnr
    INTEGER :: error
    ! ... interface
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM(modname)//'::breakinit_ks'

    !...input variables
    !   gt    : mass distribution FUNCTION
    !   xt_mg : mass of bin in mg
    !   jmax  : number of bins
    !...local variables

    DOUBLE PRECISION :: xl_d(nkr), dropradii_d(nkr), vr1_d(nkr)
    INTEGER :: ie,je,ke
    INTEGER :: i,j,k

    !...declarations for init
    INTEGER :: ip,kp,jp,kq,jq

    CHARACTER*256 file_p, file_q

    xl_d = xl_r

    ie = jbreak
    je = jbreak
    ke = jbreak

    IF (nkr == 43) file_p = 'SBM_input_43/'//'coeff_p43.dat'
    IF (nkr == 43) file_q = 'SBM_input_43/'//'coeff_q43.dat'
    IF (nkr == 33) file_p = 'SBM_input_33/'//'coeff_p_new_33.dat' ! new version 33 (taken from 43bins)
    IF (nkr == 33) file_q = 'SBM_input_33/'//'coeff_q_new_33.dat' ! new version 33   (taken from 43 bins)

    IF (my_process_is_stdio()) THEN
      unitnr=find_next_free_unit(10,999)
      WRITE(txt, '(a,i2)') 'breakinit_ks : opening coeff_p dat file'
      CALL message(modname,TRIM(txt))
      OPEN(unitnr, file=TRIM(file_p), status='old', form='formatted', iostat=error)
      IF (error /= 0) THEN
        WRITE (txt,*) 'breakinit_ks : coeff_p dat file not found'
        CALL message(modname,TRIM(txt))
        CALL finish(TRIM(modname),'error in breakinit_ks')
      END IF
      DO k=1,ke
        DO i=1,ie
          DO j=1,i
            READ (unitnr,*) kp,ip,jp,pkij(kp,ip,jp) ! pkij=[g^3*cm^3/s]
          END DO
        END DO
      END DO
      CALL message(routine,"breakinit_ks: succesfull reading coeff_p dat file")
    END IF

    IF (my_process_is_stdio()) THEN
      unitnr=find_next_free_unit(10,999)
      WRITE(txt, '(a,i2)') 'breakinit_ks : opening coeff_q dat file'
      CALL message(modname,TRIM(txt))
      OPEN(unitnr, file=TRIM(file_q), status='old', form='formatted', iostat=error)
      IF (error /= 0) THEN
        WRITE (txt,*) 'breakinit_ks : coeff_q dat file not found'
        CALL message(modname,TRIM(txt))
        CALL finish(TRIM(modname),'error in breakinit_ks')
      END IF
      DO k=1,ke
        DO j=1,je
          READ (unitnr,*) kq,jq,qkj(kq,jq)
        END DO
      END DO
      CALL message(routine,"breakinit_ks: succesfull reading coeff_q dat file")
    END IF

    dropradii_d = dropradii
    vr1_d = vr1
    DO j=1,nkr
      DO i=1,nkr
        ecoalmassm(i,j)=ecoalmass(xl_d(i), xl_d(j), dropradii_d, vr1_d, nkr)
      END DO
    END DO
    ! ... correction of coalescence efficiencies for drop collision kernels
    DO j=25,31
      ecoalmassm(nkr,j)=0.1d-29
    END DO

    RETURN

    WRITE( txt , '(a,i4)' )                                               &
    'module_fast_sbm: error opening hujisbm_data on unit, model stop' 
    CALL finish(TRIM(modname),txt)

  END SUBROUTINE breakinit_ks

  !coalescence efficiency as FUNCTION of masses
  !----------------------------------------------------------------------------+
  DOUBLE PRECISION FUNCTION ecoalmass(x1, x2, dropradii, vr1_breakup, nkr)

    IMPLICIT NONE
    INTEGER,INTENT(IN) :: nkr
    REAL(KIND=wp),INTENT(IN) :: dropradii(nkr), vr1_breakup(nkr), x1, x2

    REAL(KIND=wp) :: rho, pi, akpi, deta, dksi
    rho=1.0d0             ! [rho]=g/cm^3
    pi=3.1415927d0
    akpi=6.0d0/pi

    deta = (akpi*x1/rho)**(1.0d0/3.0d0)
    dksi = (akpi*x2/rho)**(1.0d0/3.0d0)

    ecoalmass = ecoaldiam(deta, dksi, dropradii, vr1_breakup, nkr)

    RETURN
  END FUNCTION ecoalmass
  !coalescence efficiency as FUNCTION of diameters
  !---------------------------------------------------------------------------+
  DOUBLE PRECISION FUNCTION ecoaldiam(deta,dksi,dropradii,vr1_breakup,nkr)

    IMPLICIT NONE
    INTEGER,INTENT(IN) :: nkr
    REAL(KIND=wp),INTENT(IN) :: dropradii(nkr), vr1_breakup(nkr),deta,dksi

    REAL(KIND=wp) :: dgr, dkl, rgr, rkl, q, qmin, qmax, e, x, e1, e2, sin1, cos1
    REAL(KIND=wp), PARAMETER :: one=1.0d0,eps=1.0d-30,pi=3.1415927d0

    dgr=dmax1(deta,dksi)
    dkl=dmin1(deta,dksi)

    rgr=0.5d0*dgr
    rkl=0.5d0*dkl

    q=0.5d0*(rkl+rgr)

    qmin=250.0d-4
    qmax=500.0d-4


    IF (dkl<100.0d-4) THEN
      e=1.0d0
    ELSE IF (q<qmin) THEN
      e = ecoalochs(dgr,dkl,dropradii, vr1_breakup, nkr)
    ELSE IF(q>=qmin.and.q<qmax) THEN
      x=(q-qmin)/(qmax-qmin)
      sin1=dsin(pi/2.0d0*x)
      cos1=dcos(pi/2.0d0*x)
      e1=ecoalochs(dgr, dkl, dropradii, vr1_breakup, nkr)
      e2=ecoallowlist(dgr, dkl, dropradii, vr1_breakup, nkr)
      e=cos1**2*e1+sin1**2*e2
    ELSE IF(q>=qmax) THEN
      e=ecoallowlist(dgr, dkl, dropradii, vr1_breakup, nkr)
    ELSE
      e=0.999d0
    END IF
    ecoaldiam=dmax1(dmin1(one,e),eps)

    RETURN
  END FUNCTION ecoaldiam

  !coalescence efficiency (low & list)
  !----------------------------------------------------------------------------+
  DOUBLE PRECISION FUNCTION ecoallowlist(dgr,dkl,dropradii,vr1_breakup,nkr)

    IMPLICIT NONE
    INTEGER,INTENT(IN) :: nkr
    REAL(KIND=wp),INTENT(IN) :: dropradii(nkr), vr1_breakup(nkr)
    REAL(KIND=wp),INTENT(INOUT) :: dgr, dkl

    REAL(KIND=wp) :: sigma, aka, akb, dstsc, st, sc, et, cke, qq0, qq1, qq2, ecl, w1, w2, dc
    REAL(KIND=wp), PARAMETER :: epsi=1.d-20

    ! 1 j = 10^7 g cm^2/s^2

    sigma=72.8d0    ! surface tension,[sigma]=g/s^2 (7.28e-2 n/m)
    aka=0.778d0      ! empirical constant
    akb=2.61d-4      ! empirical constant,[b]=2.61e6 m^2/j^2

    CALL collenergy(dgr,dkl,cke,st,sc,w1,w2,dc,dropradii,vr1_breakup,nkr)

    dstsc=st-sc         ! diff. of surf. energies   [dstsc] = g*cm^2/s^2
    et=cke+dstsc        ! coal. energy,             [et]    =     "

    IF (et<50.0d0) THEN    ! et < 5 uj (= 50 g*cm^2/s^2)
      qq0=1.0d0+(dkl/dgr)
      qq1=aka/qq0**2
      qq2=akb*sigma*(et**2)/(sc+epsi)
      ecl=qq1*dexp(-qq2)
    ELSE
      ecl=0.0d0
    END IF

    ecoallowlist=ecl

    RETURN
  END FUNCTION ecoallowlist

  !coalescence efficiency (beard and ochs)
  !---------------------------------------------------------------------------+
  DOUBLE PRECISION FUNCTION ecoalochs(d_l,d_s,dropradii, vr1_breakup,nkr)

    IMPLICIT NONE

    INTEGER,INTENT(IN) :: nkr
    REAL(KIND=wp),INTENT(IN) :: dropradii(nkr), vr1_breakup(nkr), d_l, d_s

    REAL(KIND=wp) :: pi, sigma, r_s, r_l, p, vtl, vts, dv, weber_number, pa1, pa2, pa3, g, x, e
    REAL(KIND=wp), PARAMETER :: fpmin=1.d-30

    pi=3.1415927d0
    sigma=72.8d0       ! surface tension [sigma] = g/s^2 (7.28e-2 n/m)
                    ! alles in cgs (1 j = 10^7 g cm^2/s^2)
    r_s=0.5d0*d_s
    r_l=0.5d0*d_l
    p=r_s/r_l

    vtl=vtbeard(d_l,dropradii, vr1_breakup,nkr)
    vts=vtbeard(d_s,dropradii, vr1_breakup,nkr)

    dv=dabs(vtl-vts)

    IF (dv<fpmin) dv=fpmin

    weber_number=r_s*dv**2/sigma

    pa1=1.0d0+p
    pa2=1.0d0+p**2
    pa3=1.0d0+p**3

    g=2**(3.0d0/2.0d0)/(6.0d0*pi)*p**4*pa1/(pa2*pa3)
    x=weber_number**(0.5d0)*g

    e=0.767d0-10.14d0*x

    ecoalochs=e

    RETURN
  END FUNCTION ecoalochs
  !ecoalochs

  !calculating the collision energy
  !------------------------------------------------------------------------------+
  SUBROUTINE collenergy(dgr,dkl,cke,st,sc,w1,w2,dc,dropradii,vr1_breakup,nkr)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: nkr
    REAL(KIND=wp),INTENT(IN) :: dropradii(nkr), vr1_breakup(nkr)
    REAL(KIND=wp),INTENT(INOUT) :: dgr, dkl, cke, st, sc, w1, w2, dc

    REAL(KIND=wp) :: pi, rho, sigma, ak10, dgka2, dgka3, v1, v2, dv, dgkb3
    REAL(KIND=wp), PARAMETER :: epsf = 1.d-30, fpmin = 1.d-30

    !external vtbeard

    pi=3.1415927d0
    rho=1.0d0            ! water density,[rho]=g/cm^3
    sigma=72.8d0         ! surf. tension,(h2o,20C)=7.28d-2 n/m
                         ! [sigma]=g/s^2
    ak10=rho*pi/12.0d0

    dgr=dmax1(dgr,epsf)
    dkl=dmax1(dkl,epsf)

    dgka2=(dgr**2)+(dkl**2)

    dgka3=(dgr**3)+(dkl**3)

    IF (dgr/=dkl) THEN
      v1=vtbeard(dgr,dropradii, vr1_breakup,nkr)
      v2=vtbeard(dkl,dropradii, vr1_breakup,nkr)
      dv=(v1-v2)
      IF (dv<fpmin) dv=fpmin
      dv=dv**2
      IF (dv<fpmin) dv=fpmin
      dgkb3=(dgr**3)*(dkl**3)
      cke=ak10*dv*dgkb3/dgka3         ! collision energy [cke]=g*cm^2/s^2
    ELSE
      cke = 0.0d0
    END IF

    st=pi*sigma*dgka2                 ! surf.energy (parent drop)
    sc=pi*sigma*dgka3**(2.0d0/3.0d0)  ! surf.energy (coal.system)

    w1=cke/(sc+epsf)                  ! weber number 1
    w2=cke/(st+epsf)                  ! weber number 2

    dc=dgka3**(1.0d0/3.0d0)           ! diam. of coal. system

    RETURN
  END SUBROUTINE collenergy
  !collenergy

  !calculating terminal velocity (beard-formula)
  !------------------------------------------------------------------------+
  DOUBLE PRECISION FUNCTION vtbeard(diam,dropradii, vr1_breakup, nkr)

    IMPLICIT NONE

    INTEGER,INTENT(IN) :: nkr
    REAL(KIND=wp),INTENT(IN) :: dropradii(nkr), vr1_breakup(nkr), diam

    INTEGER :: kr
    REAL(KIND=wp) :: aa

    aa   = diam/2.0d0           ! radius in cm

    IF (aa <= dropradii(1)) vtbeard=vr1_breakup(1)
    IF (aa > dropradii(nkr)) vtbeard=vr1_breakup(nkr)

    DO kr=1,nkr-1
      IF (aa>dropradii(kr).and.aa<=dropradii(kr+1)) THEN
        vtbeard=vr1_breakup(kr+1)
      END IF
    END DO

    RETURN
  END FUNCTION vtbeard

  SUBROUTINE ccn_init_sbm(chem_new,fccnr_con,fccnr_mar,fccnr_obs, &
                        & xland, rhocgs, zcgs )
    IMPLICIT NONE
    REAL(KIND=wp), DIMENSION(:),  INTENT(INOUT) :: chem_new
    REAL(KIND=wp), DIMENSION(:),     INTENT(IN) :: fccnr_con, &
                                              & fccnr_mar, &
                                              & fccnr_obs
    REAL(KIND=wp),                   INTENT(IN) :: xland, rhocgs, zcgs
    REAL(KIND=wp) :: factz
    INTEGER       :: kr,krr

    IF (zcgs .LE. zmin)THEN
      factz = 1.0
    ELSE
      factz=exp(-(zcgs-zmin)/z0in) !check that the height dependence is good like in mo_nwp_phy_init.f90 !
    END IF

    IF (ilognormal_modes_aerosol == 1)THEN
      ! ... generic ccn
      krr = 0
      DO kr = 1,33
        krr = krr + 1
        IF (xland > 0.9_wp) THEN !seems that my wk82 case is over sea
        !IF (xland < 0.9_wp) THEN ! pt to get continental ccn for my wk82 experiment
          chem_new(kr)=fccnr_con(krr)*factz*tune_sbmccn
          !chem_new(kr) = (fccnr_con(krr)/rhoair_max)*rhocgs ! ... distributed vertically as [#/g]
        ELSE
          chem_new(kr)=fccnr_mar(krr)*factz*tune_sbmccn
          !chem_new(kr) = (fccnr_mar(krr)/rhoair_max)*rhocgs ! ... distributed vertically as [#/g]
        END IF
        !continental anyway:
        chem_new(kr)=fccnr_con(krr)*factz*tune_sbmccn
      END DO
    ELSE
      ! ... ccn input from observation
      krr = 0
      DO kr = 1,33
        krr = krr + 1
        chem_new(kr) = fccnr_obs(krr)*factz
      END DO
    END IF
    DO kr = 1,33
      chem_new(kr)=chem_new(kr)/(rhocgs/1000.0) ! local chem_new for ccn (used in sbm): [#/cm^3]. chem_new for ccn advected by the model is [#/kg]
    END DO
  END SUBROUTINE ccn_init_sbm             



!DR  SUBROUTINE sbm_init(p_patch, p_prog_now, ext_data, p_metrics)
  SUBROUTINE sbm_init(p_patch, p_prog_now, fr_land, ddqz_z_full)

    TYPE(t_patch),    INTENT(in)    :: p_patch
    TYPE(t_nh_prog),  INTENT(inout) :: p_prog_now         ! the prognostic variables
    REAL(wp),         INTENT(IN)    :: fr_land(:,:)       ! land fraction
    REAL(wp),         INTENT(IN)    :: ddqz_z_full(:,:,:) ! cell thickness [m]

    INTEGER :: nlev
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_startidx, i_endidx    !! slices
    INTEGER :: jg, jb, jc, jk1, jk
    INTEGER :: nshift                  !< shift with respect to global grid
    REAL(KIND=wp) :: fccnr_con(iqb_s), fccnr_mar(iqb_s), fccnr_obs(iqb_s) !33?
    REAL(KIND=wp) :: dz8w, zcgs, z_full

    rl_start = 1 ! Initialization should be done for all points
    rl_end   = min_rlcell

    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

    jg = p_patch%id

    ! number of vertical levels
    nlev   = p_patch%nlev

    nshift = p_patch%nshift_total

    IF (.NOT. atm_phy_nwp_config(jg)%lsbm_warm_full) THEN !piggybacking
      !WRITE (txt,'(1X,A,i2,i2)') 'mo_sbm_util i2mom_solver,ccn_type=', &
      !atm_phy_nwp_config(jg)%cfg_2mom%i2mom_solver,atm_phy_nwp_config(jg)%cfg_2mom%ccn_type
      !CALL message('mo_sbm_util',TRIM(txt))
      !WRITE (message_text,'(1X,A,D10.3)') 'Ncn0=',ccn_coeffs%Ncn0
      !CALL message('pavel Ncn0',TRIM(message_text))
      IF (jg == 1) CALL two_moment_mcrph_init(igscp=4, msg_level=msg_level, cfg_2mom=atm_phy_nwp_config(jg)%cfg_2mom)
    END IF

    CALL warm_hucminit(fccnr_con,fccnr_mar,fccnr_obs)

    DO jb = i_startblk, i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                          &  i_startidx, i_endidx, rl_start, rl_end)
      DO jc=i_startidx,i_endidx
        z_full=0._wp
        DO jk=nlev,1,-1
          jk1 = jk + nshift !nshift is 0 in WK case
          dz8w=ddqz_z_full(jc,jk1,jb)
          zcgs=z_full+0.5*dz8w*100.0_wp
          z_full=z_full+dz8w*100.0_wp

          CALL ccn_init_sbm( &
                         & chem_new=p_prog_now%tracer(jc,jk,jb,7+iqb_s+1:7+2*iqb_s), & !jk1?
                         & fccnr_con=fccnr_con, &
                         & fccnr_mar=fccnr_mar, &
                         & fccnr_obs=fccnr_obs, &
                         & xland=fr_land(jc,jb),&
!                        & rhoair_max=0.001_wp*maxval(p_prog_now%rho(jc,:,jb)), & !is it the bottom???
                         & rhocgs=0.001_wp*p_prog_now%rho(jc,jk1,jb),   &
                         & zcgs=zcgs)
        END DO
      END DO
    END DO
  END SUBROUTINE sbm_init

END MODULE mo_sbm_util
