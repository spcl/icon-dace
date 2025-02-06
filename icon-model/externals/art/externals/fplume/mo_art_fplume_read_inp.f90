!
! mo_art_fplume_read_inp
! This module reads input files for FPLUME
! Original code by A. Folch, G. Macedonio, A. Costa (2016)
!
!
!
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

MODULE mo_art_fplume_read_inp
  USE mo_kind,                          ONLY: wp
  USE mo_art_fplume_utilities
  USE mo_impl_constants,                ONLY: MAX_CHAR_LENGTH 
  USE mo_exception,                     ONLY: finish, message, message_text 
  USE mo_art_fplume_types,              ONLY: t_fplume_phases
  USE mtime,                            ONLY: datetime, max_datetime_str_len,timedelta
  USE mo_art_data,                      ONLY: t_art_data

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: reagrn, reainp
  !
  CHARACTER(LEN=MAX_CHAR_LENGTH)                   :: mess
  CHARACTER(LEN=MAX_CHAR_LENGTH), DIMENSION(100)   :: cvoid
  INTEGER                                          :: istat,ndt0,ivoid(1),idt
  REAL(wp)                                         :: rvoid(2)
  LOGICAL                                          :: set_water_fraction     =.FALSE.
  LOGICAL                                          :: set_gas_water_fraction =.FALSE.
 
CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE reagrn(luinpname,lutgsdname, fplume_init) !  Reads granulometry from file .tgsd and creates the file .grn,                 
                  !  (accounting for aggregation if necessary)                                     
  IMPLICIT NONE                                                                                    
  !                                                                                                
  CHARACTER(LEN=*), INTENT(in) :: luinpname 
  CHARACTER(LEN=*), INTENT(in) :: lutgsdname 
  TYPE(t_fplume_phases), INTENT(inout):: fplume_init

  CHARACTER(LEN=MAX_CHAR_LENGTH) :: mess,cvoid
  INTEGER               :: istat,ic                                                                
  INTEGER(wp)           :: iaggr = 0                ! index of the first aggregating class
  REAL(wp)              :: phi_aggr,rho_aggr,frac_aggr                                             
  REAL(wp)              :: work(2),fca                                                             
  CHARACTER(LEN=MAX_CHAR_LENGTH) :: thisroutine ="mo_art_plume_read_inp:reagrn: "
  REAL(wp), ALLOCATABLE :: fip   (:)                ! fip (nc)
  REAL(wp), ALLOCATABLE :: vlimit(:)                ! vlimit(nc)
  !                                                                                                
  !*** First determine the aggregation model (if any)                                              
  !                                                                                                
  CALL get_input_cha(luinpname,'AGGREGATION','AGGREGATION_MODEL',cvoid,1,istat,mess)               
  !                                                                                                
  IF (istat/=0) THEN                               ! AGGREGATION BLOCK not found                   
     fplume_init%plume_aggregation = .FALSE.                                                                   
     fplume_init%plume_type_aggr   = 'NONE'                                                                    
  ELSE IF (TRIM(cvoid)=='NONE') THEN                                                               
     fplume_init%plume_aggregation = .FALSE.                                                                   
     fplume_init%plume_type_aggr   = 'NONE'                                                                    
  ELSE IF (TRIM(cvoid)=='CORNELL') THEN                                                            
     fplume_init%plume_aggregation = .TRUE.                                                                    
     fplume_init%plume_type_aggr   = 'CORNELL'                                                                 
  ELSE IF (TRIM(cvoid)=='PERCENTAGE') THEN                                                         
     fplume_init%plume_aggregation = .TRUE.                                                                    
     fplume_init%plume_type_aggr   = 'PERCENTAGE'                                                              
  ELSE IF (TRIM(cvoid)=='COSTA') THEN                                                              
     fplume_init%plume_aggregation = .TRUE.                                                                    
     fplume_init%plume_type_aggr   = 'COSTA'                                                                   
  ELSE                                                                                             
     CALL finish(thisroutine,'Incorrect aggregation model.')                                                    
  END IF                                                                                           
  !                                                                                                
  !*** If necessary, reads the aggregation block                                                   
  !                                                                                                
  IF (fplume_init%plume_aggregation) THEN                                                                      
     !                                                                                             
     CALL get_input_rea(luinpname,'AGGREGATION','FI_AGGREGATES',work,1,istat,mess)                 
     IF (istat>0) THEN
       WRITE(message_text,*) mess
       CALL message('',message_text)     
     ENDIF                                                           
     IF (istat<0) CALL finish(thisroutine,mess)                                                                
     phi_aggr  = work(1)                                                                           
     fplume_init%diam_aggr = 2.0_wp**(-phi_aggr)  ! diameter in mm                                             
     fplume_init%diam_aggr = 1.e-3_wp*fplume_init%diam_aggr       ! diameter in m                                              
     !                                                                                             
     CALL get_input_rea(luinpname,'AGGREGATION','DENSITY_AGGREGATES',work,1,istat,mess)            
     IF (istat>0) THEN
       WRITE(message_text,*) mess
       CALL message('',message_text)
     ENDIF                                                                
     IF (istat<0) CALL finish(thisroutine,mess)                                                                
     rho_aggr = work(1)                                                                            
     !                                                                                             
     IF (fplume_init%plume_type_aggr=='PERCENTAGE') THEN                                                       
        CALL get_input_rea(luinpname,'AGGREGATION','PERCENTAGE_(%)',work,1,istat,mess)             
        IF (istat>0) THEN
          WRITE(message_text,*) mess
          CALL message('',message_text)
        ENDIF
        IF (istat<0) CALL finish(thisroutine,mess)                                                             
        frac_aggr = 1.e-2_wp*work(1)                 ! convert to [0,1]                                
     END IF                                                                                        
     !       
     CALL get_input_rea(luinpname,'AGGREGATION','VSET_FACTOR',work,1,istat,mess)                   
     IF (istat==0) THEN                                                                            
        fplume_init%vset_aggr = work(1)                                                                        
     ELSE                                                                                          
        fplume_init%vset_aggr = 1.0_wp                                                                         
     END IF                                                                                        
     !                                                                                             
     IF (fplume_init%plume_type_aggr=='COSTA') THEN                                                            
        CALL get_input_rea(luinpname,'AGGREGATION','FRACTAL_EXPONENT',work,1,istat,mess)           
        IF (istat/=0) THEN                                                                         
           fplume_init%Dfo = 3.0_wp                                                                            
        ELSE                                                                                       
           fplume_init%Dfo = work(1)                                                                           
        END IF                                                                                     
     END IF                                                                                        
     !                                                                                             
  END IF                                                                                           
  !                                                                                                
  !*** Get the number of particles in the TGSD file                                                
  !                                                                                                
  CALL get_granulometry_nclass(lutgsdname,fplume_init%nc,istat,mess)                                           
  IF (istat>0) THEN
    WRITE(message_text,*) mess
    CALL message('',message_text)     
  ENDIF                                                              
  IF (istat<0) CALL finish(thisroutine,mess)                                                                   
  !                                                                                                
  !*** Determinates the total number of particles and allocates memory                             
  !                                                                                                
  IF (fplume_init%plume_aggregation) fplume_init%nc = fplume_init%nc + 1                                                               
  !                                                                                                
  IF (.NOT. ALLOCATED(fip)) ALLOCATE(fip   (fplume_init%nc))                                                                             
  IF (.NOT. ALLOCATED(vlimit)) ALLOCATE(vlimit(fplume_init%nc))                                                                             
  fip    = 0.0_wp                                                                                  
  vlimit = 0.0_wp                                                                                  
  IF (.NOT. ALLOCATED(fplume_init%fc)) ALLOCATE(fplume_init%fc(fplume_init%nc))
  IF (.NOT. ALLOCATED(fplume_init%rhop)) ALLOCATE(fplume_init%rhop (fplume_init%nc))
  IF (.NOT. ALLOCATED(fplume_init%diam)) ALLOCATE(fplume_init%diam(fplume_init%nc))
  IF (.NOT. ALLOCATED(fplume_init%sphe)) ALLOCATE(fplume_init%sphe(fplume_init%nc))
  IF (.NOT. ALLOCATED(fplume_init%psi)) ALLOCATE(fplume_init%psi(fplume_init%nc))
  fplume_init%fc     = 0.0_wp
  fplume_init%rhop   = 0.0_wp
  fplume_init%diam   = 0.0_wp
  fplume_init%sphe   = 0.0_wp
  fplume_init%psi    = 0.0_wp
  !                                                                                                
  !*** Reads particle properties                                                                   
  !                                                                                                
  CALL get_granulometry_value(lutgsdname,'DIAMETER',fplume_init%diam,istat,mess)                               
  IF (istat>0) THEN
    WRITE(message_text,*) mess
    CALL message('',message_text)     
  ENDIF                                                              
  IF (istat<0) CALL finish(thisroutine,mess)                                                                   
  fip (1:fplume_init%nc) = -LOG(fplume_init%diam(1:fplume_init%nc))/LOG(2.0_wp)    ! diam is in mm                                     
  fplume_init%diam(1:fplume_init%nc) = fplume_init%diam(1:fplume_init%nc)/1.e3_wp      ! convert diam to m
                                                                                                   
  CALL get_granulometry_value(lutgsdname,'DENSITY',fplume_init%rhop,istat,mess)                                
  IF (istat>0) THEN
    WRITE(message_text,*) mess
    CALL message('',message_text)
  ENDIF
  IF (istat<0) CALL finish(thisroutine,mess)                                                                   
  CALL get_granulometry_value(lutgsdname,'SPHERICITY',fplume_init%sphe,istat,mess)                             
  IF (istat>0) THEN
    WRITE(message_text,*) mess
    CALL message('',message_text)
  ENDIF
  IF (istat<0) CALL finish(thisroutine,mess)                                                                   
  CALL get_granulometry_value(lutgsdname,'FRACTION',fplume_init%fc,istat,mess)                                 
  IF (istat>0) THEN 
    WRITE(message_text,*) mess
    CALL message('',message_text)
  ENDIF
  IF (istat<0) CALL finish(thisroutine,mess)                                                                   
  !                                                                                                
  !*** Aggregated class properties                                                                 
  !                                                                                                
  IF (fplume_init%plume_aggregation) THEN                                                                      
     fplume_init%diam(fplume_init%nc) = fplume_init%diam_aggr                                                                          
     fip             (fplume_init%nc) = -LOG(1.e3_wp*fplume_init%diam_aggr)/LOG(2.0_wp)   
     fplume_init%rhop(fplume_init%nc) = rho_aggr                                                                           
     fplume_init%sphe(fplume_init%nc) = 1.0_wp                                                                             
     fplume_init%fc  (fplume_init%nc) = 0.0_wp                                                                             
  END IF                                                                                           
  !                                                                                                
  !                                                                                                
  !*** Computes the apriori aggregation (except Costa model)                                       
  !                                                                                                
  SELECT CASE(fplume_init%plume_type_aggr)                                                                     
     !                                                                                             
  CASE('NONE')                                                                                     
     !                                                                                             
     CONTINUE                                                                                      
     !                                                                                             
  CASE('PERCENTAGE')                                                                               
     !                                                                                             
     !   Computes aggregation according to a percentage                                            
     !   All classes below diam_aggr are reduced with                                              
     !   a fixed user-defined percentage                                                           
     !                                                                                             
     fca = 0.0_wp            ! fraction of mass in the aggregate class                             
     DO ic = 1,fplume_init%nc-1                                                                                
        IF (fplume_init%diam_aggr>fplume_init%diam(ic)) THEN                                                               
           fca    = fca + frac_aggr*fplume_init%fc(ic)                                                         
           fplume_init%fc(ic) = (1.0_wp-frac_aggr)*fplume_init%fc(ic)                                                      
           !                                                                                       
           IF (iaggr==0) iaggr = ic   ! index of the first aggregating class                       
           !                                                                                       
        END IF                                                                                     
     END DO                                                                                        
     fplume_init%fc(fplume_init%nc) = fca                                                                                  
                              
  CASE('CORNELL')                                                                                  
     !                                                                                             
     !   Computes aggregation according to the Cornell model                                       
     !   The aggregate class is made of                                                            
     !       50% of particles with 4<phi<=3                                                        
     !       75% of particles with 4<=phi<5                                                        
     !       90% of particles with phi>5                                                           
     !                                                                                             
     DO ic = 1,fplume_init%nc-1                                                                                
        IF (fplume_init%diam(ic)<=0.000125_wp .AND. fplume_init%diam(ic)>0.0000625_wp) THEN      
! 4<phi<3                                                                                          
           fplume_init%fc(fplume_init%nc) = fplume_init%fc(fplume_init%nc) &
                                        & + 0.5_wp*fplume_init%fc(ic)
           fplume_init%fc(ic) = 0.5_wp*fplume_init%fc(ic)                                                                  
           !                                                                                       
           IF (iaggr==0) iaggr = ic   ! index of the first aggregating class                       
           !                                                                                       
        ELSE IF (fplume_init%diam(ic)<=0.0000625_wp .AND. fplume_init%diam(ic)>=0.00003125_wp) THEN
! 4<phi<5                                                                                          
           fplume_init%fc(fplume_init%nc) = fplume_init%fc(fplume_init%nc) &
                                        & + 0.75_wp*fplume_init%fc(ic) 
           fplume_init%fc(ic) = 0.25_wp*fplume_init%fc(ic)                                                                 
           !                                                                                       
           IF (iaggr==0) iaggr = ic   ! index of the first aggregating class                       
           !                                                                                       
        ELSE IF (fplume_init%diam(ic)<0.00003125_wp) THEN    ! modified Cornell                                
           fplume_init%fc(fplume_init%nc) = fplume_init%fc(fplume_init%nc) &
                                        & + 0.9_wp*fplume_init%fc(ic)
           fplume_init%fc(ic) = 0.1_wp*fplume_init%fc(ic)                                                                  
           !                                                                                       
           IF (iaggr==0) iaggr = ic   ! index of the first aggregating class                       
           !                                                                                       
        END IF                                                                                     
     END DO                                                                                        
     !                                                                                             
  CASE('COSTA')                                                                                    
     !                                                                                             
     !*** Aggregation is computed later because it is coupled with the plume                       
     !dynamics. The grn file is also written later                                                 
     !                                                                                             
     DO ic = 1,fplume_init%nc-1                                                                                
        IF (fplume_init%diam_aggr>fplume_init%diam(ic)) THEN                                                               
           IF (iaggr==0) iaggr = ic    ! index of the first aggregating class                      
        END IF                                                                                     
     END DO                                                                                        
     IF (iaggr == 0) iaggr = fplume_init%nc-1  ! ensure at least 1 aggregating class                           
     !                                                                                             
     CONTINUE                                                                                      
     !                                                                                             
  END SELECT                                                                                       
  !                                                                                                
  RETURN                                                                                           
END SUBROUTINE reagrn 
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE reainp(luinpname,iphase, tc_exp_refdate,tc_dt_model,fplume_init)
  !********************************************************
  !*
  !*    Gets properties from the input file
  !*
  !********************************************************
  CHARACTER(LEN=*), INTENT(in)   :: luinpname
  CHARACTER(LEN=MAX_CHAR_LENGTH) :: thisroutine ="mo_art_fplume_read_inp:reainp: "      
  INTEGER, INTENT(in)            :: iphase
  TYPE(timedelta), POINTER, INTENT(in) :: &
    &  tc_dt_model                    !< Model timestep
  TYPE(datetime), POINTER, INTENT(in)  :: &
    &  tc_exp_refdate                 !< Experiment reference date configuration
  TYPE(t_fplume_phases),INTENT(inout) :: &
    &  fplume_init                     !< ART data container
  INTEGER                        :: nphases,plume_ndt
  CHARACTER(LEN=MAX_CHAR_LENGTH), DIMENSION(fplume_init%nphases+1) ::  &
    &  phase_steps
  CHARACTER(LEN=max_datetime_str_len) :: plume_beg, plume_end,phase_beg,phase_end
  INTEGER, PARAMETER :: plume_ndt_max = 100   ! maximum number of eruption phases
  CHARACTER(LEN=19)  :: plume_idt(plume_ndt_max)   ! starting time (in s) of each phase
  REAL(wp)  :: phase_Mdt,plume_Mdt(plume_ndt_max)   ! MER (kg/s) of each eruptionphase
  REAL(wp)  :: phase_Hdt,plume_Hdt(plume_ndt_max)   ! Height (m agl) of each eruption phase
  REAL(wp)  :: phase_udt,plume_udt(plume_ndt_max)   ! Exit velocity
  REAL(wp)  :: phase_Tdt,plume_Tdt(plume_ndt_max)   ! Exit temperature
  REAL(wp)  :: phase_Tvdt,plume_Tvdt(plume_ndt_max) ! Exit vapour temperature
  REAL(wp)  :: phase_Tldt,plume_Tldt(plume_ndt_max) ! Exit liquid water temperature
  REAL(wp)  :: phase_Tsdt,plume_Tsdt(plume_ndt_max) ! Exit solid water temperature
  REAL(wp)  :: phase_wvdt = 0.0_wp,plume_wvdt(plume_ndt_max)=0.0_wp  ! Exit water vapour fraction
  REAL(wp)  :: phase_wldt = 0.0_wp,plume_wldt(plume_ndt_max)=0.0_wp  ! Exit Water liquid fraction
  REAL(wp)  :: phase_wsdt = 0.0_wp,plume_wsdt(plume_ndt_max)=0.0_wp  ! Exit water solid  fraction
  REAL(wp)  :: phase_SO2,plume_SO2(plume_ndt_max)   ! MER of SO2
  REAL(wp)  :: phase_fine_ash_fraction=0.0_wp
  REAL(wp)  :: fine_ash_fraction(plume_ndt_max)=0.0_wp, ashfrac_sum
 
  !***  Gets the number of eruption intervals ndt
  !
  nphases = fplume_init%nphases

  CALL get_input_npar(luinpname,'SOURCE','EXIT_VELOCITY_(MS)',plume_ndt,istat,mess)
  IF (istat>0) THEN
    WRITE(message_text,*) mess
    CALL message('',message_text)
  ENDIF
  IF (istat<0) CALL finish(thisroutine,mess)
  !
  !***  TIME BLOCK
  !
  CALL get_input_cha(luinpname,'TIME_UTC','ERUPTION_START',cvoid,1,istat,mess)
  IF (istat>0) THEN
    WRITE(message_text,*) mess
    CALL message('',message_text)
  ENDIF
  IF (istat<0) CALL finish(thisroutine,mess)
  plume_beg = cvoid(1)
  !
  CALL get_input_cha(luinpname,'TIME_UTC','ERUPTION_END',cvoid,1,istat,mess)
  IF (istat>0) THEN
    WRITE(message_text,*) mess
    CALL message('',message_text)
  ENDIF
  IF (istat<0) CALL finish(thisroutine,mess)
  plume_end = cvoid(1)
  !
  IF (plume_ndt>1) THEN
    CALL get_input_cha(luinpname,'TIME_UTC','PHASES',cvoid,min(plume_ndt,100),istat,mess)
    IF (istat>0) THEN
      WRITE(message_text,*) mess
      CALL message('',message_text)
    ENDIF
    IF (istat<0) CALL finish(thisroutine,mess)
    DO idt = 1,plume_ndt-1
      plume_idt(idt) = cvoid(idt)
    END DO
  ENDIF
  ! construct beginning and end of plume phase
  phase_steps(1) = plume_beg
  phase_steps(nphases+1) = plume_end
  IF (nphases>1) phase_steps(2:nphases)=plume_idt(1:nphases-1)

  phase_beg = phase_steps(iphase)
  phase_end = phase_steps(iphase+1)
  !
  !***  Reads SOURCE block
  !
  CALL get_input_rea(luinpname,'SOURCE','LON_VENT',rvoid,1,istat,mess)
  IF (istat>0) THEN
    WRITE(message_text,*) mess
    CALL message('',message_text)
  ENDIF
  IF(istat<0) CALL finish(thisroutine,mess)
  fplume_init%lon=rvoid(1)
  !
  CALL get_input_rea(luinpname,'SOURCE','LAT_VENT',rvoid,1,istat,mess)
  IF (istat>0) THEN
    WRITE(message_text,*) mess
    CALL message('',message_text)
  ENDIF
  IF (istat<0) CALL finish(thisroutine,mess)
  fplume_init%lat=rvoid(1)
  !
  CALL get_input_rea(luinpname,'SOURCE','MIN_HEIGHT_FPLUME_(M)',rvoid,1,istat,mess)
  IF (istat>0) THEN
    WRITE(message_text,*) mess
    CALL message('',message_text)
  ENDIF
  IF (istat<0) CALL finish(thisroutine,mess)
  fplume_init%fplume_min_height=rvoid(1)
  !
  CALL get_input_rea(luinpname,'SOURCE','VENT_HEIGHT_(M)',rvoid,1,istat,mess)
  IF (istat>0) THEN
    WRITE(message_text,*) mess
    CALL message('',message_text)
  ENDIF
  IF (istat<0) CALL finish(thisroutine,mess)
  fplume_init%plume_zv=rvoid(1)
  !
  CALL get_input_cha(luinpname,'SOURCE','SOLVE_PLUME_FOR',fplume_init%plume_solve_for, &
                    & 1,istat,mess)
  IF(TRIM(fplume_init%plume_solve_for)=='MFR'.OR.  &
    & TRIM(fplume_init%plume_solve_for)=='mfr') THEN
     CALL get_input_rea(luinpname,'SOURCE','MFR_SEARCH_RANGE',fplume_init%plume_n_MFR,  &
                    & 2,istat,mess)
     IF (istat>0) THEN
       WRITE(message_text,*) mess
       CALL message('',message_text)
     ENDIF
     IF (istat<0) CALL finish(thisroutine,mess)
  ELSE IF (TRIM(fplume_init%plume_solve_for)=='HEIGHT' &
          & .OR. TRIM(fplume_init%plume_solve_for)=='height') THEN
     fplume_init%plume_n_MFR(1:2) = 0.0_wp
  ELSE
     CALL finish(thisroutine,'Specify the kind of plume solving strategy')
  END IF
  !
  CALL get_input_npar(luinpname,'SOURCE','HEIGHT_ABOVE_VENT_(M)',ndt0,istat,mess)
  IF (istat>0) THEN
    WRITE(message_text,*) mess
    CALL message('',message_text)
  ENDIF
  IF (istat<0) CALL finish(thisroutine,mess)
  CALL get_input_rea(luinpname,'SOURCE','HEIGHT_ABOVE_VENT_(M)',plume_Hdt,ndt0,istat,mess)
  IF (istat>0) THEN 
    WRITE(message_text,*) mess
    CALL message('',message_text)
  ENDIF
  IF (istat<0) CALL finish(thisroutine,mess)
  phase_Hdt = plume_Hdt(iphase)
  !
  CALL get_input_npar(luinpname,'SOURCE','MASS_FLOW_RATE_(KGS)',ndt0,istat,mess)
  IF (istat>0) THEN
    WRITE(message_text,*) mess
    CALL message('',message_text)
  ENDIF
  IF (istat<0) CALL finish(thisroutine,mess)
  CALL get_input_rea(luinpname,'SOURCE','MASS_FLOW_RATE_(KGS)',plume_Mdt,ndt0,istat,mess)
  IF (istat>0) THEN
    WRITE(message_text,*) mess
    CALL message('',message_text)
  ENDIF
  IF (istat<0) CALL finish(thisroutine,mess)
  phase_Mdt = plume_Mdt(iphase)
  !
  CALL get_input_npar(luinpname,'SOURCE','EXIT_VELOCITY_(MS)',ndt0,istat,mess)
  IF (istat>0) THEN
    WRITE(message_text,*) mess
    CALL message('',message_text)
  ENDIF
  IF (istat<0) CALL finish(thisroutine,mess)
  CALL get_input_rea(luinpname,'SOURCE','EXIT_VELOCITY_(MS)',plume_udt,ndt0,istat,mess)
  IF (istat>0) THEN
    WRITE(message_text,*) mess
    CALL message('',message_text)
  ENDIF
  IF (istat<0) CALL finish(thisroutine,mess)
  phase_udt = plume_udt(iphase)
  !
  CALL get_input_npar(luinpname,'SOURCE','EXIT_TEMPERATURE_(K)',ndt0,istat,mess)
  IF (istat>0) THEN
    WRITE(message_text,*) mess
    CALL message('',message_text)
  ENDIF
  IF (istat<0) CALL finish(thisroutine,mess)
  CALL get_input_rea(luinpname,'SOURCE','EXIT_TEMPERATURE_(K)',plume_Tdt,ndt0,istat,mess)
  IF (istat>0) THEN
    WRITE(message_text,*) mess
    CALL message('',message_text)
  ENDIF
  IF (istat<0) CALL finish(thisroutine,mess)
  phase_Tdt = plume_Tdt(iphase)
  !
  CALL get_input_npar(luinpname,'SOURCE','EXIT_GAS_WATER_TEMPERATURE_(K)',ndt0,istat,mess)
  ! optional
  IF (istat==0) THEN
     CALL get_input_rea(luinpname,'SOURCE','EXIT_GAS_WATER_TEMPERATURE_(K)',plume_Tvdt,ndt0,istat,mess)
     IF (istat>0) THEN 
       WRITE(message_text,*) mess
       CALL message('',message_text)
     ENDIF
     IF (istat<0) CALL finish(thisroutine,mess)
     phase_Tvdt = plume_Tvdt(iphase)
  ELSE
     phase_Tvdt = plume_Tdt(iphase)  ! Copy default values
  END IF
  !
  CALL get_input_npar(luinpname,'SOURCE','EXIT_LIQUID_WATER_TEMPERATURE_(K)',ndt0,istat,mess)
  IF (istat==0) THEN
     CALL get_input_rea(luinpname,'SOURCE','EXIT_LIQUID_WATER_TEMPERATURE_(K)',plume_Tldt,ndt0,istat,mess)
     IF (istat>0) THEN
       WRITE(message_text,*) mess
       CALL message('',message_text)
     ENDIF
     IF (istat<0) CALL finish(thisroutine,mess)
     phase_Tldt = plume_Tldt(iphase)
  ELSE
     phase_Tldt = plume_Tdt(iphase)  ! Copy default values
  END IF
  !
  CALL get_input_npar(luinpname,'SOURCE','EXIT_SOLID_WATER_TEMPERATURE_(K)',ndt0,istat,mess)
  IF (istat==0) THEN
     CALL get_input_rea(luinpname,'SOURCE','EXIT_SOLID_WATER_TEMPERATURE_(K)',plume_Tsdt,ndt0,istat,mess)
     IF (istat>0) THEN
       WRITE(message_text,*) mess
       CALL message('',message_text)
     ENDIF
     IF (istat<0) CALL finish(thisroutine,mess)
     phase_Tsdt = plume_Tsdt(iphase)
  ELSE
     phase_Tsdt = plume_Tdt(iphase)  ! Copy default values
  END IF
  !
  CALL get_input_npar(luinpname,'SOURCE','EXIT_WATER_FRACTION_(%)',ndt0,istat,mess) ! synonym of EXIT_GAS_WATER_FRACTION
  IF (istat==0) THEN
     set_water_fraction = .true.
     IF (set_gas_water_fraction) &
        CALL finish(thisroutine,'Cannot specify both EXIT_WATER_FRACTION_(%) and EXIT_GAS_WATER_FRACTION_(%)')
     CALL get_input_rea(luinpname,'SOURCE','EXIT_WATER_FRACTION_(%)',plume_wvdt,ndt0,istat,mess)
     IF (istat>0) THEN
       WRITE(message_text,*) mess
       CALL message('',message_text)
     ENDIF
     IF (istat<0) CALL finish(thisroutine,mess)
     phase_wvdt = plume_wvdt(iphase)/1.e2_wp  !Convert from %
  END IF
  !
  CALL get_input_npar(luinpname,'SOURCE','EXIT_GAS_WATER_FRACTION_(%)',ndt0,istat,mess)
  IF (istat==0) THEN
     set_gas_water_fraction = .true.
     IF (set_water_fraction) &
        CALL finish(thisroutine,                                                                  &
             &      'Cannot specify both EXIT_WATER_FRACTION_(%) and EXIT_GAS_WATER_FRACTION_(%)')
     CALL get_input_rea(luinpname,'SOURCE','EXIT_GAS_WATER_FRACTION_(%)',plume_wvdt,ndt0,istat,mess)
     IF (istat>0) THEN
       WRITE(message_text,*) mess
       CALL message('',message_text)
     ENDIF
     IF (istat<0) CALL finish(thisroutine,mess)
     phase_wvdt = plume_wvdt(iphase)/1.e2_wp   ! Convert from %
  END IF
  !
  CALL get_input_npar(luinpname,'SOURCE','EXIT_LIQUID_WATER_FRACTION_(%)',ndt0,istat,mess)
  IF (istat==0) THEN
     CALL get_input_rea(luinpname,'SOURCE','EXIT_LIQUID_WATER_FRACTION_(%)',plume_wldt,ndt0,istat,mess)
     IF (istat>0) THEN 
       WRITE(message_text,*) mess
       CALL message('',message_text)
     ENDIF
     IF (istat<0) CALL finish(thisroutine,mess)
     phase_wldt = plume_wldt(iphase)/1.e2_wp !Convert from %
  END IF
  !
  CALL get_input_npar(luinpname,'SOURCE','EXIT_SOLID_WATER_FRACTION_(%)',ndt0,istat,mess)
  IF (istat==0) THEN
     CALL get_input_rea(luinpname,'SOURCE','EXIT_SOLID_WATER_FRACTION_(%)',plume_wsdt,ndt0,istat,mess)
     IF (istat>0) THEN
       WRITE(message_text,*) mess
       CALL message('',message_text)
     ENDIF
     IF (istat<0) CALL finish(thisroutine,mess)
     phase_wsdt = plume_wsdt(iphase)/1.e2_wp
  END IF
  !
  !*** Check if at least EXIT_WATER_FRACTION or EXIT_GAS_WATER_FRACTION are set
  !
  IF (.NOT.(set_water_fraction.OR.set_gas_water_fraction)) &
       CALL finish(thisroutine,'Set either EXIT_WATER_FRACTION_(%) or EXIT_GAS_WATER_FRACTION_(%)')
  !
  CALL get_input_cha(luinpname,'SOURCE','TERMINAL_VELOCITY_MODEL',cvoid,1,istat,mess)
  IF (istat>0) THEN
    WRITE(message_text,*) mess
    CALL message('',message_text)
  ENDIF
  IF (istat<0) CALL finish(thisroutine,mess)
  IF (TRIM(cvoid(1))=='ARASTOOPOUR') THEN
     fplume_init%plume_modv = 1
  ELSE IF (TRIM(cvoid(1))=='GANSER') THEN
     fplume_init%plume_modv = 2
  ELSE IF (TRIM(cvoid(1))=='WILSON') THEN
     fplume_init%plume_modv = 3
  ELSE IF (TRIM(cvoid(1))=='DELLINO') THEN
     fplume_init%plume_modv = 4
  ELSE
     fplume_init%plume_modv = 2
     WRITE(message_text,*) 'Invalid settling velocity model. Ganser model assumed by default'
     CALL message(thisroutine,message_text)
  END IF
  !
  CALL get_input_npar(luinpname,'SOURCE','FINE_ASH_FRACTION',ndt0,istat,mess)
  IF (istat>0) THEN
    WRITE(message_text,*) mess
    CALL message('',message_text)
  ENDIF
  IF (istat<0) CALL finish(thisroutine,mess)
  IF (istat==0) THEN
    CALL get_input_rea(luinpname,'SOURCE','FINE_ASH_FRACTION',fine_ash_fraction,ndt0,istat,mess)
    IF (istat>0) THEN
      WRITE(message_text,*) mess
      CALL message('',message_text)
    ENDIF
    IF (istat<0) CALL finish(thisroutine,mess)
    phase_fine_ash_fraction = fine_ash_fraction(iphase)
  ENDIF
  !
  CALL get_input_rea(luinpname,'SOURCE','CFACTOR_ASHA',rvoid,1,istat,mess)
  IF (istat>0) THEN
    WRITE(message_text,*) mess
    CALL message('',message_text)
  ENDIF
  IF (istat<0) CALL finish(thisroutine,mess)
  fplume_init%cfactor_asha = rvoid(1)
  !
  CALL get_input_rea(luinpname,'SOURCE','CFACTOR_ASHB',rvoid,1,istat,mess)
  IF (istat>0) THEN
    WRITE(message_text,*) mess
    CALL message('',message_text)
  ENDIF
  IF (istat<0) CALL finish(thisroutine,mess)
  fplume_init%cfactor_ashb  = rvoid(1)
  !
  CALL get_input_rea(luinpname,'SOURCE','CFACTOR_ASHC',rvoid,1,istat,mess)
  IF (istat>0) THEN
    WRITE(message_text,*) mess
    CALL message('',message_text)
  ENDIF
  IF (istat<0) CALL finish(thisroutine,mess)
  fplume_init%cfactor_ashc  = rvoid(1)
  ashfrac_sum = fplume_init%cfactor_asha + fplume_init%cfactor_ashb + fplume_init%cfactor_ashc
  IF (ABS(ashfrac_sum-1.0_wp)>1e-2) THEN
    WRITE (message_text,*) ashfrac_sum
    CALL message(thisroutine,'sum of ash mode fractions is not equal to 1: '//message_text)
  ENDIF
  !
  CALL get_input_npar(luinpname,'SOURCE','MER_SO2',ndt0,istat,mess)
  IF (istat>0) THEN
    WRITE(message_text,*) mess
    CALL message('',message_text)
  ENDIF
  IF (istat<0) CALL finish(thisroutine,mess)
  CALL get_input_rea(luinpname,'SOURCE','MER_SO2',plume_SO2,ndt0,istat,mess)
  IF (istat>0) THEN
    WRITE(message_text,*) mess
    CALL message('',message_text)
  ENDIF
  IF (istat<0) CALL finish(thisroutine,mess)
  phase_SO2 = plume_SO2(iphase)
  !
  !*** Reads CONSTANTS block
  !
  CALL get_input_rea(luinpname,'CONSTANTS','Cw',rvoid,1,istat,mess)
  IF (istat==0) fplume_init%Cw = rvoid(1)
  !
  CALL get_input_rea(luinpname,'CONSTANTS','Cp',rvoid,1,istat,mess)
  IF (istat==0) fplume_init%Cp = rvoid(1)
  !
  CALL get_input_rea(luinpname,'CONSTANTS','Ca',rvoid,1,istat,mess)
  IF (istat==0) fplume_init%Ca = rvoid(1)
  !
  !*** Reads the MODELLING block (optional, all variables have default values if
  !not found)
  !
  CALL get_input_cha(luinpname,'MODELLING','WIND_COUPLING',cvoid,1,istat,mess)
  IF (istat==0) THEN
     IF (TRIM(cvoid(1))=='NO'.OR.TRIM(cvoid(1))=='no') THEN
       fplume_init%plume_wind_coupling = .false.
     ENDIF
  END IF
  !
  CALL get_input_cha(luinpname,'MODELLING','AIR_MOISTURE',cvoid,1,istat,mess)
  IF (istat==0) THEN
     IF (TRIM(cvoid(1))=='NO'.OR.TRIM(cvoid(1))=='no') THEN
       fplume_init%plume_moist_air = .false.
     ENDIF
  END IF
  !
  CALL get_input_cha(luinpname,'MODELLING','LATENT_HEAT',cvoid,1,istat,mess)
  IF (istat==0) THEN
     IF (TRIM(cvoid(1))=='NO'.OR.TRIM(cvoid(1))=='no') THEN 
       fplume_init%plume_latent_heat = .false.
     ENDIF
  END IF
  !
  CALL get_input_cha(luinpname,'MODELLING','REENTRAINMENT',cvoid,1,istat,mess)
  IF (istat==0) THEN
     IF (TRIM(cvoid(1))=='NO'.OR.TRIM(cvoid(1))=='no') THEN
       fplume_init%plume_reentrainment = .false.
     ENDIF
  end if
  !
  CALL get_input_rea (luinpname,'MODELLING','xi',rvoid,1,istat,mess)
  IF (istat==0) fplume_init%plume_xi = rvoid(1)
  !
  CALL get_input_rea (luinpname,'MODELLING','zmin_wind',rvoid,1,istat,mess)
  IF (istat==0) fplume_init%plume_zmin_wind = rvoid(1)
  !
  CALL get_input_rea (luinpname,'MODELLING','c_umbrella',rvoid,1,istat,mess)
  IF (istat==0) fplume_init%plume_c_umbrella = rvoid(1)
  !
  CALL get_input_cha(luinpname,'MODELLING','a_s',cvoid,1,istat,mess)
  IF (istat==0) THEN
     IF (TRIM(cvoid(1))=='CONSTANT'.OR.TRIM(cvoid(1))=='constant') THEN
        fplume_init%plume_type_as = 'CONSTANT'
        CALL get_input_rea(luinpname,'MODELLING','a_s',rvoid,2,istat,mess)
        IF (istat==0) THEN
           fplume_init%plume_a_s_jet   = rvoid(1)
           fplume_init%plume_a_s_plume = rvoid(2)
        ELSE
           CALL finish(thisroutine,mess)
        END IF
     ELSE IF (TRIM(cvoid(1))=='KAMINSKI-R'.OR.TRIM(cvoid(1))=='kaminski-r') THEN
        fplume_init%plume_type_as = 'KAMINSKI-R'
     ELSE IF (TRIM(cvoid(1))=='KAMINSKI-C'.OR.TRIM(cvoid(1))=='kaminski-c') THEN
        fplume_init%plume_type_as = 'KAMINSKI-C'
     ELSE IF (TRIM(cvoid(1))=='OLD'.OR.TRIM(cvoid(1))=='old') THEN
        fplume_init%plume_type_as = 'OLD'
     ELSE
        CONTINUE
     END IF
  END IF
  !
  CALL get_input_cha(luinpname,'MODELLING','a_v',cvoid,1,istat,mess)
  IF (istat==0) THEN
     IF (TRIM(cvoid(1))=='CONSTANT'.OR.TRIM(cvoid(1))=='constant') THEN
        fplume_init%plume_type_av = 'CONSTANT'
        CALL get_input_rea(luinpname,'MODELLING','a_v',rvoid,1,istat,mess)
        IF (istat==0) THEN
           fplume_init%plume_a_v = rvoid(1)
        ELSE
           CALL finish(thisroutine,mess)
        END IF
     ELSE IF (TRIM(cvoid(1))=='TATE'.OR.TRIM(cvoid(1))=='tate') THEN
        fplume_init%plume_type_av = 'TATE'
     ELSE
        CONTINUE
     END IF
  END IF

  CALL setpsi(fplume_init)   ! Calculates psi(nc) depending on the model
  !
  IF (iphase==1) ALLOCATE(fplume_init%phase(nphases))
  CALL fplume_init%phase(iphase)%init(phase_beg,phase_end,phase_SO2,                  &
                            & phase_Hdt, phase_Mdt, phase_udt, phase_Tdt,phase_Tvdt,phase_Tldt,  &
                            & phase_Tsdt, phase_wvdt,phase_wldt,phase_wsdt,tc_exp_refdate,       &
                            & tc_dt_model,phase_fine_ash_fraction)

  RETURN
END SUBROUTINE reainp
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE setpsi(fplume_init)
  !*********************************************************************
  !*
  !*      Calculates the particle shape factor psi depending on the
  !*      velocity model
  !*
  !*      modv = 1   ARASTOPOUR. psi = 1 (not used)
  !*      modv = 2   GANSER      psi = sphericity
  !*      modv = 3   WILSON      psi = (b+c)/2a    a>b>c semi-axes
  !*      modv = 4   DELLINO     psi = sphericity/circularity
  !*
  !**********************************************************************
  IMPLICIT NONE
  TYPE(t_fplume_phases), INTENT(inout):: fplume_init
  !
  INTEGER    :: ic
  REAL(wp)   :: pi,gama,circula
  !
  !***  Initializations
  pi = 4.0_wp*atan(1.0_wp)
  !
  !***  Computes psi
  IF (fplume_init%plume_modv==1) THEN          ! Arastopour
     DO ic=1,fplume_init%nc
        fplume_init%psi(ic) = 1.0_wp
     ENDDO
  ELSE If (fplume_init%plume_modv==2) THEN     ! Ganser
     DO ic=1,fplume_init%nc
        fplume_init%psi(ic) = fplume_init%sphe(ic)
     ENDDO
  ELSE IF (fplume_init%plume_modv==3) THEN     ! Wilson
     DO ic = 1,fplume_init%nc
        CALL get_gama(gama,fplume_init%diam(ic),fplume_init%sphe(ic))  ! Get a/c
        IF(gama>=1.0_wp) THEN               ! oblate
           fplume_init%psi(ic) = 0.5_wp*(1.0_wp+1.0_wp/gama)
        ELSE                                ! prolate
           fplume_init%psi(ic) = gama
        ENDIF
     ENDDO
  ELSE IF (fplume_init%plume_modv==4) THEN    ! Dellino
     DO ic = 1,fplume_init%nc
        circula = 1.69_wp
        fplume_init%psi(ic) = fplume_init%sphe(ic)/circula
     ENDDO
  ENDIF
  RETURN
END SUBROUTINE setpsi
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE get_gama(gama,diam,sphe)
  !**************************************************************************
  !*
  !*     Gets gama = a/c
  !*
  !*     NOTE: In all cases it is assumed that particles fall as
  !*     prolate ellipsoids
  !*
  !*            a = b < c  prolate   (gama < 1)
  !*
  !*     The inversion of the area of the ellipsoid is done numerically.
  !*     Area given by:
  !*
  !*     A = 2*pi*(a**2 + c**2*e/tan(e) )   e = acos(gama)
  !*     d = 2*c*gama**(2/3)               (prolate)
  !*
  !*     NOTE: particle diameter is multiplied by a factor. It does not affect
  !*           results (a/c) and is done to facilitate convergence and prevent
  !*           propagation of rounding errors (e.g. for micron size particles
  !*           diam of the order 1d-6 rised to 2 or 3)
  !*
  !***************************************************************************
  IMPLICIT NONE
  REAL(wp) :: gama,diam,sphe
  !
  INTEGER    :: iiter,niter
  REAL(wp)   :: d,pi,gmin,gmax,Ao,toler,e
  REAL(wp)   :: Vp,Ap
  CHARACTER(LEN=MAX_CHAR_LENGTH) :: thisroutine ="mo_art_fplume_read_inp:get_gama: "  
  !
  !***   Initializations
  d     = diam*1.e4_wp         ! see NOTE
  niter = 1000
  toler = 1.e-6_wp
  gmin  = 1.e-3_wp
  gmax  = 1.0_wp
  !
  !***   Volume and area
  pi = 4.0_wp*ATAN(1.0_wp)
  Vp = 4.0_wp*pi*((0.5_wp*d)**3)/3.0_wp
  Ap = (pi**(1.0_wp/3.0_wp))*((6.0_wp*Vp)**(2.0_wp/3.0_wp))/sphe
  !
  !***   Iterates
  DO iiter = 1,niter
     gama = 0.5_wp*(gmin+gmax)
     e    = ACOS(gama)
     Ao   = 0.5_wp*pi*d*d*(gama**(-4.0_wp/3.0_wp))*(gama*gama + (e/TAN(e)))
     IF (Ao<Ap) THEN
        gmax = gama
     ELSE
        gmin = gama
     ENDIF
     IF((iiter>1).AND.(ABS(Ao-Ap).lt.toler)) EXIT 
     IF (iiter==niter) THEN
       WRITE(message_text,*) 'convergence not achieved'
       CALL message(thisroutine,message_text)
       EXIT
     ENDIF
  ENDDO
  RETURN 
END SUBROUTINE get_gama
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_fplume_read_inp
