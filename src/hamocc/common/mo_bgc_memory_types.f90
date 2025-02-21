! Definition of variables and allocation of memory
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


MODULE mo_bgc_memory_types

  USE mo_kind, ONLY   : wp
  USE mo_exception, ONLY      : message
  
  
  USE mo_param1_bgc, ONLY: n_bgctra, natm, npowtra, nbgctend,nbgcflux, &
    & nsedtra, nsed_diag

  USE mo_control_bgc, ONLY: rmasko, bgc_nproma, bgc_zlevs
  USE mo_hamocc_nml, ONLY: l_cpl_co2, & 
    & isac,ks,ksp,dzs,porwat
 
  USE mo_memory_agg, ONLY: naggdiag
  
  IMPLICIT NONE

  PUBLIC
 
  !-------------------------------------------------------------------------
  TYPE t_bgc_memory

    INTEGER,  POINTER  :: kbo(:)   !< k-index of bottom layer (2d)

    REAL(wp), POINTER :: bgctra(:,:,:)
    REAL(wp), POINTER :: bgctend(:,:,:)
    REAL(wp), POINTER :: bgcflux(:,:)
    REAL(wp), POINTER :: solco2(:)

    REAL(wp), POINTER :: co3 (:,:)
    REAL(wp), POINTER :: hi  (:,:)
    REAL(wp), POINTER :: aksp(:,:)

    REAL(wp), POINTER :: swr_frac (:,:)
    REAL(wp), POINTER :: meanswr(:,:)
    REAL(wp), POINTER :: akw3(:,:)
    REAL(wp), POINTER :: akb3(:,:)
    REAL(wp), POINTER :: ak13(:,:)
    REAL(wp), POINTER :: ak23(:,:)
    REAL(wp), POINTER :: aks3(:,:)
    REAL(wp), POINTER :: akf3(:,:)
    REAL(wp), POINTER :: ak1p3(:,:)
    REAL(wp), POINTER :: ak2p3(:,:)
    REAL(wp), POINTER :: ak3p3(:,:)
    REAL(wp), POINTER :: aksi3(:,:)
    REAL(wp), POINTER :: satoxy(:,:)
    REAL(wp), POINTER :: satn2(:)
    REAL(wp), POINTER :: satn2o(:)
    REAL(wp), POINTER :: atdifv(:)       ! not used
    REAL(wp), POINTER :: suppco2(:)      ! not used
    REAL(wp), POINTER :: sedfluxo(:,:)
    REAL(wp), POINTER :: aksurf(:,:)    !> dissociation constants for DIC,bor and water 
                                              !> at the sea surface, no pressure dependency 
    REAL(wp), POINTER :: atm(:,:)  

    REAL(wp), POINTER :: co2flux(:)     !> sea-air C-flux, ingetrated over whole simulation period,  not used
    REAL(wp), POINTER :: o2flux(:)      !> sea-air O2-flux, ingetrated over whole simulation period, not used
    REAL(wp), POINTER :: n2flux(:)      !> sea-air N2-flux, ingetrated over whole simulation period, not used
    REAL(wp), POINTER :: n2oflux(:)     !> sea-air N2-flux, ingetrated over whole simulation period, not used
    REAL(wp), POINTER :: co2trans(:)    !> transfer coefficient for CO2 atmosph/ocean, not used
    REAL(wp), POINTER :: co2conc(:)     ! not used
    REAL(wp), POINTER :: co2flux_cpl(:) ! not used

    REAL(wp), POINTER :: bolay(:) !<  height of bottom cell
    REAL(wp), POINTER :: expoor(:)  ! not used
    REAL(wp), POINTER :: expoca(:)  ! not used
    REAL(wp), POINTER :: exposi(:)  ! not used
    REAL(wp), POINTER :: strahl(:)

    REAL(wp), POINTER :: alar1max(:)  ! not used
    REAL(wp), POINTER :: TSFmax(:)    ! not used
    REAL(wp), POINTER :: TMFmax(:)    ! not used
    
    REAL(wp), POINTER :: wpoc(:,:)        ! depth-dependent detritus settling speed, global variable
    REAL(wp), POINTER :: wdust(:,:)     ! depth-dependent dust settling speed
    REAL(wp), POINTER :: wopal(:,:)     ! daily sinking speed of opal (namelist parameter)
    REAL(wp), POINTER :: wcal(:,:)     ! daily sinking speed of cal (namelist parameter)
    REAL(wp), POINTER :: ws_agg(:,:)    ! daily sinking speed of aggregates (namelist parameter)
    

  END TYPE t_bgc_memory
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  TYPE t_sediment_memory
  
    REAL(wp), POINTER :: sedlay (:,:,:)
    REAL(wp), POINTER :: powtra (:,:,:)
    REAL(wp), POINTER :: sedtend(:,:,:)

    REAL(wp), POINTER :: sedhpl(:,:)
    REAL(wp), POINTER :: seddenit(:)

    REAL(wp), POINTER :: silpro(:)
    REAL(wp), POINTER :: prorca(:)
    REAL(wp), POINTER :: prcaca(:)
    REAL(wp), POINTER :: produs(:)

    REAL(wp), POINTER :: burial(:,:)

    ! pown2bud closes the mass balance for the alkalinity for biogenic induced changes in N2 in sediment
    REAL(wp), POINTER :: pown2bud(:,:)
    ! powh2obud closes the mass balance for oxygen
    REAL(wp), POINTER :: powh2obud(:,:)

  END TYPE t_sediment_memory
  !-------------------------------------------------------------------------
  
  
  !-------------------------------------------------------------------------
  TYPE t_aggregates_memory
      REAL(wp),POINTER :: av_dp(:,:),              &  ! mean primary particle diameter
                       &  av_rho_p(:,:),           &  ! mean primary particle density
                       &  df_agg(:,:),             &  ! fractal dimension of aggregates
                       &  b_agg(:,:),              &  ! aggregate number distribution slope
                       &  Lmax_agg(:,:),           &  ! maximum diamater of aggregates
                       &  ws_agg(:,:),             &  ! aggregate mean sinking velocity
                       &  stickiness_agg(:,:),     &  ! mean aggregate stickiness
                       &  stickiness_frustule(:,:),&  ! frustule stickiness
                       &  dynvis(:,:),             &  ! molecular dynamic viscosity
                       &  av_rhof_V(:,:) !,          &  ! volume-weighted aggregate density
!                       &  av_d_c(:,:),             &  ! concentration-weighted mean diameter of aggs
!                       &  av_por_V(:,:)               ! volume-weighted aggregate porosity

    REAL(wp), POINTER :: aggdiag(:,:,:)    ! 3d concentration EU

  END TYPE t_aggregates_memory
  !-------------------------------------------------------------------------

    
  TYPE(t_bgc_memory), ALLOCATABLE, TARGET :: bgc_local_memory(:)
  TYPE(t_sediment_memory), ALLOCATABLE, TARGET :: sediment_local_memory(:)
  TYPE(t_aggregates_memory), ALLOCATABLE, TARGET :: aggregates_memory(:)
  
!   TYPE(t_bgc_memory), POINTER :: bgc_memory
!   TYPE(t_sediment_memory), POINTER :: sediment_memory

  INTEGER:: bgc_memory_copies = 1

CONTAINS


  !-------------------------------------------------------------------------
  SUBROUTINE allocate_all_bgc_memory
  
    INTEGER :: i
    CHARACTER(LEN=*), PARAMETER :: routine_name = 'allocate_all_bgc_memory'


    !
    ! Allocate allocate_bgc_memory_types
    !
    CALL message(TRIM(routine_name), 'allocate_bgc_memory_types')
    ALLOCATE(bgc_local_memory(0:bgc_memory_copies-1))
    DO i=0,bgc_memory_copies-1
      CALL allocate_bgc_memory_types(bgc_local_memory(i))
    ENDDO
    !
    ! Allocate memory : biology
    ! Already done in allocate_bgc_memory_types
    !
  !   CALL message(TRIM(routine_name), 'alloc_mem_biomod')
  !   CALL ALLOC_MEM_BIOMOD

     !
    ! Allocate memory : sediment
    !
    CALL message(TRIM(routine_name), 'alloc_mem_sedmnt' )
    ALLOCATE(sediment_local_memory(0:bgc_memory_copies-1))
    DO i=0,bgc_memory_copies-1
      CALL allocate_mem_sediment_types(sediment_local_memory(i))
    ENDDO
    !
    ! Allocate memory : inorganic carbon cycle
    !
    ! done in allocate_bgc_memory_types
!     CALL message(TRIM(routine_name), 'alloc_mem_carbch')
     
    CALL message(TRIM(routine_name), 'alloc_mem_aggregates')
    ALLOCATE(aggregates_memory(0:bgc_memory_copies-1))
    DO i=0,bgc_memory_copies-1
      CALL allocate_mem_aggregates_types(aggregates_memory(i))
    ENDDO
    
  END SUBROUTINE allocate_all_bgc_memory
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE allocate_mem_aggregates_types(local_aggregates_memory)
    TYPE(t_aggregates_memory) :: local_aggregates_memory
    
        !-----------------------------------------------------------------------
     !>
     !! Initialization/allocation fields
     !! Called in ini_bgc after read_namelist
     !!
     ! allocate memory space for aggregate properties
     ALLOCATE(local_aggregates_memory%av_dp(bgc_nproma,bgc_zlevs))
     ALLOCATE(local_aggregates_memory%av_rho_p(bgc_nproma,bgc_zlevs))
     ALLOCATE(local_aggregates_memory%df_agg(bgc_nproma,bgc_zlevs))
     ALLOCATE(local_aggregates_memory%b_agg(bgc_nproma,bgc_zlevs))
     ALLOCATE(local_aggregates_memory%Lmax_agg(bgc_nproma,bgc_zlevs))
     ALLOCATE(local_aggregates_memory%stickiness_agg(bgc_nproma,bgc_zlevs))
     ALLOCATE(local_aggregates_memory%stickiness_frustule(bgc_nproma,bgc_zlevs))
     ALLOCATE(local_aggregates_memory%av_rhof_V(bgc_nproma,bgc_zlevs))
!     ALLOCATE(av_d_C(bgc_nproma,bgc_zlevs))
!     ALLOCATE(av_por_V(bgc_nproma,bgc_zlevs))

     ALLOCATE(local_aggregates_memory%aggdiag(bgc_nproma,bgc_zlevs,naggdiag))

     ! mean sinking velocity
     ALLOCATE(local_aggregates_memory%ws_agg(bgc_nproma,bgc_zlevs))

     ! molecular dynamic viscosity
     ALLOCATE(local_aggregates_memory%dynvis(bgc_nproma,bgc_zlevs))

     local_aggregates_memory%av_dp = 0._wp
     local_aggregates_memory%av_rho_p = 0._wp
     local_aggregates_memory%df_agg = 0._wp
     local_aggregates_memory%b_agg = 0._wp
     local_aggregates_memory%Lmax_agg = 0._wp
     local_aggregates_memory%stickiness_agg = 0._wp
     local_aggregates_memory%stickiness_frustule = 0._wp
     local_aggregates_memory%av_rhof_V = 0._wp
!     av_d_C = 0._wp
!     av_por_V = 0._wp
     local_aggregates_memory%aggdiag = 0._wp

     !$ACC ENTER DATA COPYIN(local_aggregates_memory, local_aggregates_memory%aggdiag)

  END SUBROUTINE allocate_mem_aggregates_types
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE deallocate_mem_aggregates_types(local_aggregates_memory)
    TYPE(t_aggregates_memory) :: local_aggregates_memory

     DEALLOCATE(local_aggregates_memory%av_dp)
     DEALLOCATE(local_aggregates_memory%av_rho_p)
     DEALLOCATE(local_aggregates_memory%df_agg)
     DEALLOCATE(local_aggregates_memory%b_agg)
     DEALLOCATE(local_aggregates_memory%Lmax_agg)
     DEALLOCATE(local_aggregates_memory%stickiness_agg)
     DEALLOCATE(local_aggregates_memory%stickiness_frustule)
     DEALLOCATE(local_aggregates_memory%aggdiag)
     DEALLOCATE(local_aggregates_memory%ws_agg)
     DEALLOCATE(local_aggregates_memory%dynvis)
     DEALLOCATE(local_aggregates_memory%av_rhof_V)
!     DEALLOCATE(av_d_C)
!     DEALLOCATE(av_por_V)

  END SUBROUTINE deallocate_mem_aggregates_types
  !-------------------------------------------------------------------------
 
  !-------------------------------------------------------------------------
  SUBROUTINE allocate_bgc_memory_types(bgc_mem_instance)
    TYPE(t_bgc_memory) :: bgc_mem_instance
    
    CALL allocate_memory_carbch(bgc_mem_instance)
    CALL alloc_mem_biomod_types(bgc_mem_instance)

    !$ACC ENTER DATA COPYIN(bgc_mem_instance, bgc_mem_instance%atm, bgc_mem_instance%satn2, bgc_mem_instance%satn2o) &
    !$ACC   COPYIN(bgc_mem_instance%solco2, bgc_mem_instance%bgctra, bgc_mem_instance%hi, bgc_mem_instance%co3) &
    !$ACC   COPYIN(bgc_mem_instance%aksp, bgc_mem_instance%ak13, bgc_mem_instance%ak23, bgc_mem_instance%akb3) &
    !$ACC   COPYIN(bgc_mem_instance%akw3, bgc_mem_instance%aksi3, bgc_mem_instance%ak1p3, bgc_mem_instance%ak2p3) &
    !$ACC   COPYIN(bgc_mem_instance%ak3p3, bgc_mem_instance%aks3, bgc_mem_instance%akf3, bgc_mem_instance%satoxy) &
    !$ACC   COPYIN(bgc_mem_instance%bgctend, bgc_mem_instance%kbo, bgc_mem_instance%bolay, bgc_mem_instance%strahl) &
    !$ACC   COPYIN(bgc_mem_instance%swr_frac, bgc_mem_instance%meanswr, bgc_mem_instance%bgcflux, bgc_mem_instance%aksurf) &
    !$ACC   COPYIN(bgc_mem_instance%wpoc, bgc_mem_instance%wcal, bgc_mem_instance%wopal, bgc_mem_instance%wdust) &
    !$ACC   COPYIN(bgc_mem_instance%sedfluxo)

!     bgc_memory => bgc_local_memory(1)
  
  END SUBROUTINE allocate_bgc_memory_types
  !------------------------------------------------------------------------- 
  
  !-------------------------------------------------------------------------
  SUBROUTINE allocate_memory_carbch(bgc_mem_instance)
    
    TYPE(t_bgc_memory) :: bgc_mem_instance

    ALLOCATE (bgc_mem_instance%bgctra(bgc_nproma,bgc_zlevs,n_bgctra))
    
    ALLOCATE (bgc_mem_instance%bgctend(bgc_nproma,bgc_zlevs,nbgctend))
    bgc_mem_instance%bgctend = 0._wp

    ALLOCATE (bgc_mem_instance%bgcflux(bgc_nproma,nbgcflux))
    bgc_mem_instance%bgcflux = 0._wp

    ALLOCATE (bgc_mem_instance%hi(bgc_nproma,bgc_zlevs))

    ALLOCATE (bgc_mem_instance%co3(bgc_nproma,bgc_zlevs))

    ALLOCATE (bgc_mem_instance%solco2(bgc_nproma))
    bgc_mem_instance%solco2(:) = rmasko

    ALLOCATE (bgc_mem_instance%satn2o(bgc_nproma))
    bgc_mem_instance%satn2o(:) = rmasko

    ALLOCATE (bgc_mem_instance%satoxy(bgc_nproma,bgc_zlevs))
    bgc_mem_instance%satoxy(:,:) = rmasko

    ALLOCATE (bgc_mem_instance%satn2(bgc_nproma))
    bgc_mem_instance%satn2(:) = rmasko

    ALLOCATE (bgc_mem_instance%sedfluxo(bgc_nproma,npowtra))
    bgc_mem_instance%sedfluxo(:,:) = 0._wp

    ALLOCATE (bgc_mem_instance%aksp(bgc_nproma,bgc_zlevs))

    ALLOCATE (bgc_mem_instance%aks3(bgc_nproma,bgc_zlevs))

    ALLOCATE (bgc_mem_instance%akf3(bgc_nproma,bgc_zlevs))

    ALLOCATE (bgc_mem_instance%ak1p3(bgc_nproma,bgc_zlevs))

    ALLOCATE (bgc_mem_instance%ak2p3(bgc_nproma,bgc_zlevs))

    ALLOCATE (bgc_mem_instance%ak3p3(bgc_nproma,bgc_zlevs))

    ALLOCATE (bgc_mem_instance%aksi3(bgc_nproma,bgc_zlevs))

    ALLOCATE (bgc_mem_instance%ak23(bgc_nproma,bgc_zlevs))

    ALLOCATE (bgc_mem_instance%ak13(bgc_nproma,bgc_zlevs))

    ALLOCATE (bgc_mem_instance%akb3(bgc_nproma,bgc_zlevs))

    ALLOCATE (bgc_mem_instance%akw3(bgc_nproma,bgc_zlevs))

    ALLOCATE (bgc_mem_instance%swr_frac(bgc_nproma,bgc_zlevs))
    bgc_mem_instance%swr_frac=1._wp
    ALLOCATE (bgc_mem_instance%meanswr(bgc_nproma,bgc_zlevs))
    bgc_mem_instance%meanswr=0._wp
    ALLOCATE (bgc_mem_instance%aksurf(bgc_nproma,10))
    bgc_mem_instance%aksurf(:,:) = rmasko
    ALLOCATE (bgc_mem_instance%atm(bgc_nproma,natm))
    ALLOCATE (bgc_mem_instance%atdifv(bgc_nproma))
    ALLOCATE (bgc_mem_instance%suppco2(bgc_nproma))
    bgc_mem_instance%suppco2(:) = 0.0_wp

    ALLOCATE (bgc_mem_instance%co2flux(bgc_nproma))
    bgc_mem_instance%co2flux(:) = 0.0_wp
    ALLOCATE (bgc_mem_instance%o2flux(bgc_nproma))
    bgc_mem_instance%o2flux(:) = 0.0_wp
    ALLOCATE (bgc_mem_instance%n2flux(bgc_nproma))
    bgc_mem_instance%n2flux(:) = 0.0_wp
    ALLOCATE (bgc_mem_instance%n2oflux(bgc_nproma))
    bgc_mem_instance%n2oflux(:) = 0.0_wp

    IF (l_cpl_co2) THEN
      ALLOCATE (bgc_mem_instance%co2trans(bgc_nproma))
      bgc_mem_instance%co2trans(:) = 0.0_wp
      ALLOCATE (bgc_mem_instance%co2conc(bgc_nproma))
      bgc_mem_instance%co2conc(:) = 0.0_wp
      ALLOCATE (bgc_mem_instance%co2flux_cpl(bgc_nproma))
      bgc_mem_instance%co2flux_cpl(:) = 0.0_wp
    ENDIF

    !ToDo: probably there is abetter place for these allocations
    ALLOCATE (bgc_mem_instance%kbo(bgc_nproma))
        
    ALLOCATE (bgc_mem_instance%wpoc(bgc_nproma,bgc_zlevs))
    ALLOCATE (bgc_mem_instance%wdust(bgc_nproma,bgc_zlevs))
    ALLOCATE (bgc_mem_instance%wopal(bgc_nproma,bgc_zlevs))
    ALLOCATE (bgc_mem_instance%wcal(bgc_nproma,bgc_zlevs))
    ALLOCATE (bgc_mem_instance%ws_agg(bgc_nproma,bgc_zlevs))
        
  END SUBROUTINE allocate_memory_carbch
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE deallocate_all_bgc_memory_types()
  
    INTEGER :: i
 
    DO i=1,bgc_memory_copies-1
      CALL deallocate_bgc_local_memory_types(bgc_local_memory(i))
      CALL deallocate_sediment_local_mem(sediment_local_memory(i))
      CALL deallocate_mem_aggregates_types(aggregates_memory(i))
    ENDDO
    
    DEALLOCATE(bgc_local_memory)
    DEALLOCATE(sediment_local_memory)
    DEALLOCATE(aggregates_memory)
 
  END SUBROUTINE deallocate_all_bgc_memory_types
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE deallocate_bgc_local_memory_types(bgc_mem_instance)
    
    TYPE(t_bgc_memory) :: bgc_mem_instance
  
    DEALLOCATE (bgc_mem_instance%bgctra)
    
    DEALLOCATE (bgc_mem_instance%bgctend)
 
    DEALLOCATE (bgc_mem_instance%bgcflux)

    DEALLOCATE (bgc_mem_instance%hi)

    DEALLOCATE (bgc_mem_instance%co3)

    DEALLOCATE (bgc_mem_instance%solco2)

    DEALLOCATE (bgc_mem_instance%satn2o)

    DEALLOCATE (bgc_mem_instance%satoxy)

    DEALLOCATE (bgc_mem_instance%satn2)
 
    DEALLOCATE (bgc_mem_instance%sedfluxo)

    DEALLOCATE (bgc_mem_instance%aksp)

    DEALLOCATE (bgc_mem_instance%aks3)

    DEALLOCATE (bgc_mem_instance%akf3)

    DEALLOCATE (bgc_mem_instance%ak1p3)

    DEALLOCATE (bgc_mem_instance%ak2p3)

    DEALLOCATE (bgc_mem_instance%ak3p3)

    DEALLOCATE (bgc_mem_instance%aksi3)

    DEALLOCATE (bgc_mem_instance%ak23)

    DEALLOCATE (bgc_mem_instance%ak13)

    DEALLOCATE (bgc_mem_instance%akb3)

    DEALLOCATE (bgc_mem_instance%akw3)

    DEALLOCATE (bgc_mem_instance%swr_frac)
    DEALLOCATE (bgc_mem_instance%meanswr)
    DEALLOCATE (bgc_mem_instance%aksurf)
    DEALLOCATE (bgc_mem_instance%atm)
    DEALLOCATE (bgc_mem_instance%atdifv)
    DEALLOCATE (bgc_mem_instance%suppco2)

    DEALLOCATE (bgc_mem_instance%co2flux)
    DEALLOCATE (bgc_mem_instance%o2flux)
    DEALLOCATE (bgc_mem_instance%n2flux)
    DEALLOCATE (bgc_mem_instance%n2oflux)

    IF (l_cpl_co2) THEN
      DEALLOCATE (bgc_mem_instance%co2trans)
      DEALLOCATE (bgc_mem_instance%co2conc)
      DEALLOCATE (bgc_mem_instance%co2flux_cpl)
    ENDIF
    
    ! dealloc_mem_biomod_types
    DEALLOCATE (bgc_mem_instance%expoor)
    DEALLOCATE (bgc_mem_instance%expoca)
    DEALLOCATE (bgc_mem_instance%exposi)
    DEALLOCATE (bgc_mem_instance%kbo)
    DEALLOCATE (bgc_mem_instance%strahl)
    DEALLOCATE (bgc_mem_instance%bolay)
    DEALLOCATE (bgc_mem_instance%alar1max)
    DEALLOCATE (bgc_mem_instance%TSFmax)
    DEALLOCATE (bgc_mem_instance%TMFmax)

  END SUBROUTINE deallocate_bgc_local_memory_types
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  SUBROUTINE alloc_mem_biomod_types(bgc_mem_instance)
    USE mo_control_bgc
    USE mo_param1_bgc
    
    TYPE(t_bgc_memory) :: bgc_mem_instance

    ALLOCATE (bgc_mem_instance%expoor(bgc_nproma))
    bgc_mem_instance%expoor(:)=0._wp
    ALLOCATE (bgc_mem_instance%expoca(bgc_nproma))
    bgc_mem_instance%expoca(:)=0._wp
    ALLOCATE (bgc_mem_instance%exposi(bgc_nproma))
    bgc_mem_instance%exposi(:)=0._wp
    ALLOCATE (bgc_mem_instance%kbo(bgc_nproma))
    ALLOCATE (bgc_mem_instance%strahl(bgc_nproma))
    bgc_mem_instance%strahl=200._wp
    ALLOCATE (bgc_mem_instance%bolay(bgc_nproma))
    ALLOCATE (bgc_mem_instance%alar1max(bgc_nproma))
    ALLOCATE (bgc_mem_instance%TSFmax(bgc_nproma))
    ALLOCATE (bgc_mem_instance%TMFmax(bgc_nproma))

  END SUBROUTINE alloc_mem_biomod_types
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE allocate_mem_sediment_types(sediment_local_mem)
    TYPE(t_sediment_memory) :: sediment_local_mem

    CALL allocate_local_mem_sediment(sediment_local_mem)
!     sediment_memory => sediment_local_memory(1)
    
  END SUBROUTINE allocate_mem_sediment_types
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE allocate_local_mem_sediment(sediment_local_mem)
    TYPE(t_sediment_memory) :: sediment_local_mem
  
    ALLOCATE (sediment_local_mem%sedlay(bgc_nproma,ks,nsedtra))
    sediment_local_mem%sedlay(:,:,:) = 0._wp

    ALLOCATE (sediment_local_mem%sedhpl(bgc_nproma,ks))
    ALLOCATE (sediment_local_mem%burial(bgc_nproma,nsedtra))
    sediment_local_mem%burial(:,:) = 0._wp
    ALLOCATE (sediment_local_mem%powtra(bgc_nproma,ks,npowtra))
    ALLOCATE (sediment_local_mem%silpro(bgc_nproma))
    sediment_local_mem%silpro(:) = 0._wp
    ALLOCATE (sediment_local_mem%prorca(bgc_nproma))
    sediment_local_mem%prorca(:) = 0._wp

    ALLOCATE (sediment_local_mem%prcaca(bgc_nproma))
    sediment_local_mem%prcaca(:) = 0._wp
    ALLOCATE (sediment_local_mem%produs(bgc_nproma))
    sediment_local_mem%produs(:) = 0._wp
    ALLOCATE (sediment_local_mem%pown2bud(bgc_nproma,ks))
    sediment_local_mem%pown2bud(:,:) = 0._wp
    ALLOCATE (sediment_local_mem%powh2obud(bgc_nproma,ks))
    sediment_local_mem%powh2obud(:,:) = 0._wp
    ALLOCATE (sediment_local_mem%sedtend(bgc_nproma,ks,nsed_diag))
    sediment_local_mem%sedtend(:,:,:)=0._wp
    ALLOCATE (sediment_local_mem%seddenit(bgc_nproma))
    sediment_local_mem%seddenit(:)=0._wp

    !$ACC ENTER DATA COPYIN(sediment_local_mem, sediment_local_mem%silpro, sediment_local_mem%produs) &
    !$ACC   COPYIN(sediment_local_mem%prcaca, sediment_local_mem%prorca, sediment_local_mem%burial) &
    !$ACC   COPYIN(sediment_local_mem%sedlay, sediment_local_mem%powtra, sediment_local_mem%sedhpl) &
    !$ACC   COPYIN(sediment_local_mem%powh2obud, sediment_local_mem%pown2bud, sediment_local_mem%sedtend)

  END SUBROUTINE allocate_local_mem_sediment
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  SUBROUTINE deallocate_sediment_local_mem(sediment_local_mem)
    TYPE(t_sediment_memory) :: sediment_local_mem
    
    DEALLOCATE (sediment_local_mem%sedlay)

    DEALLOCATE (sediment_local_mem%sedhpl)
    DEALLOCATE (sediment_local_mem%burial)
    DEALLOCATE (sediment_local_mem%powtra)
    DEALLOCATE (sediment_local_mem%silpro)
    DEALLOCATE (sediment_local_mem%prorca)

    DEALLOCATE (sediment_local_mem%prcaca)
    DEALLOCATE (sediment_local_mem%produs)
    DEALLOCATE (sediment_local_mem%pown2bud)
    DEALLOCATE (sediment_local_mem%powh2obud)
    DEALLOCATE (sediment_local_mem%sedtend)
    DEALLOCATE (sediment_local_mem%seddenit)

  END SUBROUTINE deallocate_sediment_local_mem
  !-------------------------------------------------------------------------
END MODULE mo_bgc_memory_types
  !-------------------------------------------------------------------------
