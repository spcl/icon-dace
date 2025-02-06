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

! @brief Variables for sediment modules
!
! declaration and memory allocation

MODULE mo_sedmnt

  USE mo_kind, ONLY        : wp
  USE mo_exception, ONLY      : message, finish
  USE mo_param1_bgc, ONLY  : nsedtra, npowtra, n_bgctra, nsed_diag
  USE mo_control_bgc, ONLY: dtbgc, bgc_nproma, bgc_zlevs 
  USE mo_hamocc_nml, ONLY : isac,ks,ksp,dzs,porwat
  USE mo_memory_bgc, ONLY :  sinkspeed_dust
  USE mo_bgc_constants, ONLY: g,rhoref_water
  USE mo_bgc_memory_types
  USE mo_fortran_tools, ONLY: set_acc_host_or_device

  IMPLICIT NONE

  PUBLIC

  REAL(wp), ALLOCATABLE :: seddw(:)   ! global variable
  REAL(wp), ALLOCATABLE :: porsol(:)  ! global variable
  REAL(wp), ALLOCATABLE :: porwah(:)  ! global variable
  REAL(wp), ALLOCATABLE :: pors2w(:)  ! global variable

  REAL(wp), ALLOCATABLE :: seddzi(:)  ! global variable
  REAL(wp), ALLOCATABLE :: z_sed(:)   ! global variable

  REAL(wp), ALLOCATABLE :: zcoefsu(:), zcoeflo(:)

  REAL(wp) :: sedict, calcon, rno3, o2ut, ansed, sedac, sedifti
  REAL(wp) :: calcwei, opalwei, orgwei
  REAL(wp) :: calcdens, opaldens, orgdens, claydens
  REAL(wp) :: calfa, oplfa, orgfa, clafa, solfu
  REAL(wp) :: sred_sed, disso_op,disso_cal
  REAL(wp) :: silsat, o2thresh
CONTAINS

SUBROUTINE sediment_bottom

 REAL(wp) :: sumsed, dustd1, dustd2
 INTEGER:: k
  !

  seddzi(1) = 500._wp

  DO k = 1, ks
     porsol(k) = 1._wp - porwat(k)
     IF ( k == 1 ) porwah(k) = 0.5_wp * (1._wp + porwat(1))
     IF ( k >= 2 ) porwah(k) = 0.5_wp * (porwat(k) + porwat(k-1))
     ! inverse thickness of sediment layer
     seddzi(k+1) = 1._wp / dzs(k+1)
     ! define sediment layer thickness [m]
     seddw(k) = 0.5_wp * (dzs(k) + dzs(k+1))
  END DO

  ! ******************************************************************
  ! the following section is to include the SEDIMENT ACCELERATION
  ! mechanism to accelerate the sediment:

  sedac = 1._wp / REAL(isac, wp)

  ! determine total solid sediment thickness sumsed
  ! and reduced sediment volume

  sumsed = 0._wp

  DO k = 1, ks
     porwat(k) = porwat(k) * sedac
     porsol(k) = porsol(k) * sedac
     pors2w(k) = 1._wp/porwat(k) -1._wp ! ratio of porsol/porwat
     sumsed = sumsed + seddw(k)
  ENDDO

  ! depth of sediment layers for output

  z_sed(1) = dzs(1)

  DO k = 2, ks
     z_sed(k) = z_sed(k-1) + dzs(k)
  ENDDO

  ! determine reduced diffusion rate sedict,
  ! and scaling for bottom layer ventilation, sedifti

  sedict=1.e-9_wp*dtbgc                        ! units? m/s2 ?
  sedifti = sedict / (sumsed**2)               ! not used
  sedict=sedict*sedac

  zcoefsu(0) = 0.0_wp
  DO  k = 1, ks
     ! sediment diffusion coefficient * 1/dz * fraction of pore water at half depths
     zcoefsu(k  ) = -sedict * seddzi(k) * porwah(k)
     ! the lowerdiffusive flux of layer k is identical to
     ! the upper diff. flux of layer k+1
     zcoeflo(k-1) = -sedict * seddzi(k) * porwah(k)
  END DO

  ! diffusion coefficient for bottom sediment layer
  zcoeflo(ks) = 0.0_wp

  ! ******************************************************************
  !
  ! densities etc. for SEDIMENT SHIFTING
  !
  ! define wei(ght) of calcium carbonate, opal, and poc [kg/kmol]
  !
  calcwei = 100._wp ! 40+12+3*16 kg/kmol C
  opalwei = 60._wp  ! 28 + 2*16  kg/kmol Si
  orgwei  = 30._wp  ! from 12 kg/kmol * 2.5 POC[kg]/DW[kg]
                    ! after Alldredge, 1998: POC(g)/DW(g) = 0.4 of diatom marine snow, size 1mm3
  !
  ! define dens(ities) of caco3, opal, caco3 [kg/m3]
  !
  calcdens = 2600._wp
  opaldens = 2200._wp
  orgdens  = 1000._wp
  claydens = 2600._wp !quartz
  !
  ! define volumes occupied by solid constituents [m3/kmol]
  !
  calfa = calcwei/calcdens         ! 3.85e-2
  oplfa = opalwei/opaldens         ! 2.73e-2
  orgfa = orgwei/orgdens           ! 3.0e-2
  clafa = 1._wp/claydens           ! 3.85e-4 (clay is calculated in kg/m3)
  !
  ! determine total solid column integrated sediment volume  (1-d)
  !
  solfu = 0._wp
  !
  DO k=1,ks
     solfu = solfu + seddw(k)*porsol(k)
  ENDDO

  dustd1 = 0.0001_wp !cm = 1 um, boundary between clay and silt
  dustd2 = dustd1*dustd1
  sinkspeed_dust = (g * 86400._wp / 18._wp                      &  ! g * sec per day / 18. | js: Stoke's law for small particles
         &   * (claydens - rhoref_water) / 1.567_wp * 1000._wp  &  ! excess density / dyn. visc. | -> cm/s to m/day
         &   * dustd2 * 1.e-4_wp)                              ! *diameter**2 |*1000 *1.e-4?

  !$ACC ENTER DATA COPYIN(pors2w, seddw, porsol, porwat, zcoefsu, zcoeflo)


END SUBROUTINE



SUBROUTINE  ini_bottom(local_bgc_mem, start_idx, end_idx, klevs, pddpo, lacc)

 TYPE(t_bgc_memory), POINTER :: local_bgc_mem
 REAL(wp), INTENT(IN):: pddpo(bgc_nproma,bgc_zlevs)

 INTEGER, INTENT(IN) :: start_idx, end_idx
 INTEGER,  TARGET::klevs(bgc_nproma)
 LOGICAL, INTENT(IN), OPTIONAL :: lacc

 LOGICAL :: lzacc
 INTEGER ::  j, k, kpke

  CALL set_acc_host_or_device(lzacc, lacc)

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO j = start_idx, end_idx
   k=klevs(j)
        local_bgc_mem%kbo(j) = 1
        local_bgc_mem%bolay(j) = 0._wp
        IF(k.gt.0)THEN
        IF ( pddpo(j,k) > EPSILON(0.5_wp) ) THEN
           local_bgc_mem%bolay(j) = pddpo(j,k)
           local_bgc_mem%kbo(j) = k
        ENDIF
        ENDIF
  END DO

 ! find depth of last wet layer
  ! and evaluate global min depth of last layer

  ! bolaymin is neanilges, as is only compuetd on an arbitery local memory
!   bolaymin=8000._wp
  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO j = start_idx, end_idx
   kpke=klevs(j)
   IF(kpke>0)THEN ! always valid for MPIOM
   !$ACC LOOP SEQ
   DO k = kpke-1, 1, -1
      IF ( pddpo(j,k) > EPSILON(0.5_wp) .AND. pddpo(j,k+1) <= EPSILON(0.5_wp) ) THEN
         local_bgc_mem%bolay(j) = pddpo(j,k) ! local thickness of bottom layer
         local_bgc_mem%kbo(j) = k
!          bolaymin = MIN(bolaymin,local_bgc_mem%bolay(j))
      ENDIF
    END DO
   ENDIF
  END DO
  !$ACC END PARALLEL
  
 END SUBROUTINE
 
 SUBROUTINE ALLOC_MEM_SEDMNT

    ALLOCATE (seddzi(ksp))
    ALLOCATE (z_sed(ks))
    ALLOCATE (seddw(ks))
    ALLOCATE (porsol(ks))
    ALLOCATE (porwah(ks))
    ALLOCATE (pors2w(ks))
    ALLOCATE (zcoefsu(0:ks))
    ALLOCATE (zcoeflo(0:ks))

  END SUBROUTINE ALLOC_MEM_SEDMNT

END MODULE mo_sedmnt
