! Contains the variables to set up the ocean model.
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

! DSL definitions 





























!---------------------
! block definitions




!---------------------
! mappings



















!---------------------
! connectivity









!---------------------
! generic types




!---------------------
! shortcuts







































!---------------------
! Upper-lower case

MODULE mo_hamocc_types

  USE mo_kind,                ONLY: wp, sp
  USE mo_impl_constants,      ONLY: land, land_boundary, boundary, sea_boundary, sea,  &
    & success, max_char_length, min_dolic,               &
    & full_coriolis, beta_plane_coriolis,                &
    & f_plane_coriolis, zero_coriolis, halo_levels_ceiling
  USE mo_math_types,         ONLY: t_cartesian_coordinates,      &
    & t_geographical_coordinates
  USE mo_ocean_tracer_transport_types, ONLY: t_tracer_collection


  PUBLIC :: t_hamocc_diag
  PUBLIC :: t_hamocc_monitor
  PUBLIC :: t_hamocc_state
  PUBLIC :: t_hamocc_tend
  PUBLIC :: t_hamocc_sed
  PUBLIC :: t_hamocc_bcond
  PUBLIC :: t_hamocc_prog
  PUBLIC :: t_hamocc_agg

  !
  TYPE t_onCells_Pointer_3d_wp
    REAL(wp),  POINTER, DIMENSION(:,:,:) :: p  ! pointer to 3D array
  END TYPE t_onCells_Pointer_3d_wp

  TYPE t_hamocc_monitor
    REAL(wp), POINTER :: phosy(:)
    REAL(wp), POINTER :: phosy_cya(:)
    REAL(wp), POINTER :: grazing(:)
    REAL(wp), POINTER :: omex90(:)
    REAL(wp), POINTER :: calex90(:)
    REAL(wp), POINTER :: opex90(:)
    REAL(wp), POINTER :: omex1000(:)
    REAL(wp), POINTER :: calex1000(:)
    REAL(wp), POINTER :: opex1000(:)
    REAL(wp), POINTER :: omex2000(:)
    REAL(wp), POINTER :: calex2000(:)
    REAL(wp), POINTER :: opex2000(:)
    REAL(wp), POINTER :: net_co2_flux(:)
    REAL(wp), POINTER :: sfalk(:)
    REAL(wp), POINTER :: sfdic(:)
    REAL(wp), POINTER :: sfphos(:)
    REAL(wp), POINTER :: sfsil(:)
    REAL(wp), POINTER :: sfnit(:)
    REAL(wp), POINTER :: carbinv(:)
    REAL(wp), POINTER :: bacfra(:)
    REAL(wp), POINTER :: phymor(:)
    REAL(wp), POINTER :: zoomor(:)
    REAL(wp), POINTER :: exudz(:)
    REAL(wp), POINTER :: graton(:)
    REAL(wp), POINTER :: exud(:)
    REAL(wp), POINTER :: n2fix(:)
    REAL(wp), POINTER :: wcdenit(:)
    REAL(wp), POINTER :: remina(:)
    REAL(wp), POINTER :: remins(:)
    REAL(wp), POINTER :: seddenit(:)
    REAL(wp), POINTER :: cyaldet(:)
    REAL(wp), POINTER :: cyaldoc(:)
    REAL(wp), POINTER :: delsil(:)
    REAL(wp), POINTER :: delcar(:)
    REAL(wp), POINTER :: zalkn2(:)
    ! N-cycle
    REAL(wp), POINTER :: phosy_nh4(:)
    REAL(wp), POINTER :: phosy_cya_nh4(:)
    REAL(wp), POINTER :: sfnh4(:)
    REAL(wp), POINTER :: net_nh3_flux(:)
    REAL(wp), POINTER :: wc_nitri_no2(:)
    REAL(wp), POINTER :: wc_nitri_nh4(:)
    REAL(wp), POINTER :: wc_dnrn(:)
    REAL(wp), POINTER :: wc_dnra(:)
    REAL(wp), POINTER :: wc_anammox(:)


  END TYPE t_hamocc_monitor

  TYPE t_hamocc_diag
    !--------------------------------------------
    ! dimension: (nproma,n_zlev, nblks_e)
    REAL(wp), POINTER ::  hi(:,:,:)          
    REAL(wp), POINTER ::  co3(:,:,:)     
    !--------------------------------------------
  END TYPE t_hamocc_diag
  !

  TYPE t_hamocc_sed
    !--------------------------------------------
    ! dimension: (nproma, ks, nblks_e)
    REAL(wp), POINTER ::  so12(:,:,:)       
    REAL(wp), POINTER ::  sc12(:,:,:)       
    REAL(wp), POINTER ::  ssil(:,:,:)       
    REAL(wp), POINTER ::  ster(:,:,:)       
    REAL(wp), POINTER ::  pwic(:,:,:)       
    REAL(wp), POINTER ::  pwal(:,:,:)       
    REAL(wp), POINTER ::  pwph(:,:,:)       
    REAL(wp), POINTER ::  pwox(:,:,:)       
    REAL(wp), POINTER ::  pwsi(:,:,:)       
    REAL(wp), POINTER ::  pwfe(:,:,:)       
    REAL(wp), POINTER ::  pwh2s(:,:,:)       
    REAL(wp), POINTER ::  pwn2(:,:,:)       
    REAL(wp), POINTER ::  pwno3(:,:,:)       
    REAL(wp), POINTER ::  sedhi(:,:,:)       
    REAL(wp), POINTER ::  pwh2ob(:,:,:)       
    REAL(wp), POINTER ::  pwn2b(:,:,:)       
    REAL(wp), POINTER ::  bo12(:,:)       
    REAL(wp), POINTER ::  bc12(:,:)       
    REAL(wp), POINTER ::  bsil(:,:)       
    REAL(wp), POINTER ::  bter(:,:)       

    ! N-cycle
    REAL(wp), POINTER ::  pwnh4(:,:,:)
    REAL(wp), POINTER ::  pwno2(:,:,:)   
 
    INTEGER, POINTER ::  kbo(:,:)       
    REAL(wp), POINTER ::  bolay(:,:)       
  END TYPE t_hamocc_sed

  TYPE t_hamocc_tend
    !--------------------------------------------
    ! dimension: (nproma,n_zlev, nblks_e)
    REAL(wp), POINTER ::  npp(:,:,:)           
    REAL(wp), POINTER ::  phoc(:,:,:)           
    REAL(wp), POINTER ::  cyloss(:,:,:)           
    REAL(wp), POINTER ::  graz(:,:,:)       
    REAL(wp), POINTER ::  graton(:,:,:)       
    REAL(wp), POINTER ::  exud(:,:,:)       
    REAL(wp), POINTER ::  exudz(:,:,:)       
    REAL(wp), POINTER ::  zoomor(:,:,:)       
    REAL(wp), POINTER ::  phymor(:,:,:)       
    REAL(wp), POINTER ::  nfix(:,:,:)       
    REAL(wp), POINTER ::  remina(:,:,:)       
    REAL(wp), POINTER ::  remins(:,:,:)       
    REAL(wp), POINTER ::  reminn(:,:,:)       
    REAL(wp), POINTER ::  bacfra(:,:,:)       
    REAL(wp), POINTER ::  co2mr(:,:)       
    REAL(wp), POINTER ::  pco2(:,:)       
    REAL(wp), POINTER ::  cflux(:,:)       
    REAL(wp), POINTER ::  nflux(:,:)       
    REAL(wp), POINTER ::  n2oflux(:,:)       
    REAL(wp), POINTER ::  dmsflux(:,:)       
    REAL(wp), POINTER ::  oflux(:,:)       
    REAL(wp), POINTER ::  nfixd(:,:)       
    REAL(wp), POINTER ::  delsil(:,:,:)       
    REAL(wp), POINTER ::  delcar(:,:,:)       
    REAL(wp), POINTER ::  prcaca(:,:)       
    REAL(wp), POINTER ::  produs(:,:)       
    REAL(wp), POINTER ::  prorca(:,:)       
    REAL(wp), POINTER ::  silpro(:,:)       
    REAL(wp), POINTER ::  coex90(:,:)       
    REAL(wp), POINTER ::  opex90(:,:)       
    REAL(wp), POINTER ::  calex90(:,:)       
    REAL(wp), POINTER ::  coex1000(:,:)       
    REAL(wp), POINTER ::  opex1000(:,:)       
    REAL(wp), POINTER ::  calex1000(:,:)       
    REAL(wp), POINTER ::  coex2000(:,:)       
    REAL(wp), POINTER ::  opex2000(:,:)       
    REAL(wp), POINTER ::  calex2000(:,:)       
    REAL(wp), POINTER ::  orginp(:,:)   ! for now constant value read via nml    
    REAL(wp), POINTER ::  silinp(:,:)   ! later riverine input possible     
    REAL(wp), POINTER ::  calinp(:,:)   ! via mo_bgc_bcond    
  !  REAL(wp), POINTER ::  nitinp(:,:)     ! nitrogen input 
    REAL(wp), POINTER ::  h2obudget(:,:,:)       
    REAL(wp), POINTER ::  n2budget(:,:,:)       
    REAL(wp), POINTER ::  sedflic(:,:)       
    REAL(wp), POINTER ::  sedflal(:,:)       
    REAL(wp), POINTER ::  sedflph(:,:)       
    REAL(wp), POINTER ::  sedflox(:,:)       
    REAL(wp), POINTER ::  sedflsi(:,:)       
    REAL(wp), POINTER ::  sedflfe(:,:)       
    REAL(wp), POINTER ::  sedfln2(:,:)       
    REAL(wp), POINTER ::  sedflno3(:,:)       
    REAL(wp), POINTER ::  sedflh2s(:,:)       
    REAL(wp), POINTER ::  sedro2(:,:,:)       
    REAL(wp), POINTER ::  sedrn(:,:,:)       
    REAL(wp), POINTER ::  sedrs(:,:,:)       
    REAL(wp), POINTER ::  dmsprod(:,:,:)       
    REAL(wp), POINTER ::  dmsbac(:,:,:)       
    REAL(wp), POINTER ::  dmsuv(:,:,:)       
    REAL(wp), POINTER ::  euexp(:,:,:)       
    REAL(wp), POINTER ::  plim(:,:,:)       
    REAL(wp), POINTER ::  flim(:,:,:)       
    REAL(wp), POINTER ::  nlim(:,:,:)       
    REAL(wp), POINTER ::  aksp(:,:,:)       
    REAL(wp), POINTER ::  akb(:,:,:)       
    REAL(wp), POINTER ::  akw(:,:,:)       
    REAL(wp), POINTER ::  ak1(:,:,:)       
    REAL(wp), POINTER ::  ak2(:,:,:)       
    REAL(wp), POINTER ::  aksi(:,:,:)
    REAL(wp), POINTER ::  aks(:,:,:)
    REAL(wp), POINTER ::  akf(:,:,:)
    REAL(wp), POINTER ::  ak1p(:,:,:)
    REAL(wp), POINTER ::  ak2p(:,:,:)
    REAL(wp), POINTER ::  ak3p(:,:,:)
    REAL(wp), POINTER ::  satoxy(:,:,:)       
    REAL(wp), POINTER ::  satn2(:,:)       
    REAL(wp), POINTER ::  satn2o(:,:)       
    REAL(wp), POINTER ::  solco2(:,:)       
    REAL(wp), POINTER ::  aou(:,:,:)       
    REAL(wp), POINTER ::  cTlim(:,:,:)       
    REAL(wp), POINTER ::  cLlim(:,:,:)       
    REAL(wp), POINTER ::  cPlim(:,:,:)       
    REAL(wp), POINTER ::  cFlim(:,:,:)       
    REAL(wp), POINTER ::  o2min(:,:)       
    REAL(wp), POINTER ::  zo2min(:,:)       
    REAL(wp), POINTER ::  h2sprod(:,:,:)       
    REAL(wp), POINTER ::  h2sloss(:,:,:)       
    REAL(wp), POINTER ::  lysocline(:,:)
    REAL(wp), POINTER ::  nitrogeninp(:,:)

    ! N-cycle
    REAL(wp), POINTER ::  nh3flux(:,:)
    REAL(wp), POINTER ::  gppnh4(:,:,:)
    REAL(wp), POINTER ::  cyapro(:,:,:)
    REAL(wp), POINTER ::  ammox(:,:,:)
    REAL(wp), POINTER ::  nitox(:,:,:)
    REAL(wp), POINTER ::  dnrn(:,:,:)
    REAL(wp), POINTER ::  dnra(:,:,:)
    REAL(wp), POINTER ::  anam(:,:,:)
    REAL(wp), POINTER ::  sedflnh4(:,:)       
    REAL(wp), POINTER ::  sedflno2(:,:)  
    REAL(wp), POINTER ::  sedammox(:,:,:)       
    REAL(wp), POINTER ::  sednitox(:,:,:)       
    REAL(wp), POINTER ::  sedanam(:,:,:)
    REAL(wp), POINTER ::  seddnrn(:,:,:)       
    REAL(wp), POINTER ::  seddnra(:,:,:)       
    REAL(wp), POINTER ::  sednrn2(:,:,:)     
    REAL(wp), POINTER ::  wdust(:,:,:)
    REAL(wp), POINTER ::  wpoc(:,:,:)
    REAL(wp), POINTER ::  wopal(:,:,:)
    REAL(wp), POINTER ::  wcal(:,:,:)

    TYPE(t_hamocc_monitor) :: monitor
  END TYPE t_hamocc_tend

 TYPE t_hamocc_agg  ! mm: agg
    !--------------------------------------------
    ! dimension: (nproma,n_zlev, nblks_e)
    REAL(wp), POINTER ::  wsagg(:,:,:)
    REAL(wp), POINTER ::  avdp(:,:,:)
    REAL(wp), POINTER ::  avrhop(:,:,:)
    REAL(wp), POINTER ::  sticka(:,:,:)
    REAL(wp), POINTER ::  stickf(:,:,:)
    REAL(wp), POINTER ::  lmaxagg(:,:,:)
    REAL(wp), POINTER ::  dfagg(:,:,:)
    REAL(wp), POINTER ::  avdc(:,:,:)
    REAL(wp), POINTER ::  bagg(:,:,:)
    REAL(wp), POINTER ::  dynvis(:,:,:)
    REAL(wp), POINTER ::  avrhof(:,:,:)
    REAL(wp), POINTER ::  avpor(:,:,:)
  END TYPE t_hamocc_agg

  TYPE t_hamocc_prog
    REAL(wp), POINTER :: tracer(:,:,:,:)
    TYPE(t_tracer_collection) :: tracer_collection
    TYPE(t_onCells_Pointer_3d_wp), ALLOCATABLE :: tracer_ptr(:) 
  END TYPE t_hamocc_prog

! array of states
   TYPE t_hamocc_state

    TYPE(t_hamocc_diag) :: p_diag
    TYPE(t_hamocc_tend) :: p_tend
    TYPE(t_hamocc_sed)  :: p_sed
    TYPE(t_hamocc_agg)  :: p_agg
    TYPE(t_hamocc_prog), ALLOCATABLE :: p_prog(:)
    
  END TYPE t_hamocc_state

  TYPE t_hamocc_bcond
   REAL(wp), POINTER:: dusty(:,:)       !  index1=1,nproma, nblks_e
   REAL(wp), POINTER:: nitro(:,:)       !  index1=1,nproma, nblks_e
   REAL(wp), POINTER:: prorca(:,:)       !  index1=1,nproma, nblks_e
   REAL(wp), POINTER:: prcaca(:,:)       !  index1=1,nproma, nblks_e
   REAL(wp), POINTER:: silpro(:,:)       !  index1=1,nproma, nblks_e
   REAL(wp), POINTER:: produs(:,:)       !  index1=1,nproma, nblks_e
  END TYPE t_hamocc_bcond

END MODULE 


