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

! Contains the types to set up the wave model.

MODULE mo_wave_types

  USE mo_kind,                ONLY: wp, vp
  USE mo_var_list,            ONLY: t_var_list_ptr
  USE mo_fortran_tools,       ONLY: t_ptr_2d3d, t_ptr_2d_int

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: t_wave_prog
  PUBLIC :: t_wave_source
  PUBLIC :: t_wave_diag
  PUBLIC :: t_wave_state
  PUBLIC :: t_wave_state_lists

  TYPE t_wave_prog
    REAL(wp), POINTER, CONTIGUOUS :: &
    tracer(:,:,:,:) => NULL()
    !! wave energy (spectral bins) over frequencies and directions (nproma, nlev, nblks_c, ntracer) [m^2 ?]
    TYPE(t_ptr_2d3d), ALLOCATABLE :: tracer_ptr(:) !< pointer array: one pointer for each tracer
  END TYPE t_wave_prog


  ! source function state vector object
  !
  TYPE t_wave_source
    REAL(wp), POINTER, CONTIGUOUS :: &
      &  sl(:,:,:),           & ! total source function                    (nproma,nblks_c,ntracer) (-)
      &  fl(:,:,:)              ! diagonal matrix of functional derivative (nproma,nblks_c,ntracer) (-)

    INTEGER, POINTER, CONTIGUOUS ::  &
      &  llws(:,:,:)            ! 1 - where sinput is positive (nproma,nblks_c,ntracer) (-)

    TYPE(t_ptr_2d3d),   ALLOCATABLE :: sl_ptr(:)   !< pointer array: one pointer for each tracer
    TYPE(t_ptr_2d3d),   ALLOCATABLE :: fl_ptr(:)   !< pointer array: one pointer for each tracer
    TYPE(t_ptr_2d_int), ALLOCATABLE :: llws_ptr(:) !< pointer array: one pointer for each tracer
  END TYPE t_wave_source


  ! diagnostic variables state vector
  !
  TYPE t_wave_diag
    REAL(wp), POINTER, CONTIGUOUS :: &
      &  gv_c(:,:,:),         & ! group velocity                    (nproma,nblks_c,nfreqs)  (m/s)
      &  gv_e(:,:,:),         & ! group velocity                    (nproma,nblks_e,nfreqs)  (m/s)
      &  gvn_e(:,:,:,:),      & ! orthogonal normal group velocity  (nproma,nlev,nblks_e,nfreqs)  (m/s)
      &  gvt_e(:,:,:,:),      & ! tangential group velocity         (nproma,nlev,nblks_e,nfreqs)  (m/s)
      &  alphaj(:,:),         & ! jonswap alpha                     (nproma,nblks_c)         (-)
      &  fp(:,:),             & ! jonswap peak frequency            (nproma,nblks_c)         (hz)
      &  et(:,:,:),           & ! jonswap spectra                   (nproma,nblks_c,nfreqs)  (-)
      &  flminfr(:,:,:),      & ! the minimum value in spectral bins for a given frequency (nproma,nblks_c,nfreqs)
      &  f1mean(:,:),         & ! mean frequency based on f-moment  (nproma,nblks_c)
      &  wave_num_c(:,:,:),   & ! wave number at cell centers as a function of
                                ! circular frequency and water depth (nproma,nblks_c,nfreqs) (1/m)
      &  wave_num_e(:,:,:),   & ! wave number at cell edges as a function of 
                                ! circular frequency and water depth (nproma,nblks_e,nfreqs) (1/m)
      &  akmean(:,:),         & ! mean wavenumber based on sqrt(1/k)-moment  (nproma,nblks_c) (1/m)
      &  xkmean(:,:),         & ! mean wavenumber based on sqrt(k)-moment    (nproma,nblks_c) (1/m)
      &  ustar(:,:),          & ! friction velocity                          (nproma,nblks_c) (m/s)
      &  z0(:,:),             & ! roughness length                           (nproma,nblks_c) (m)
      &  tauhf1(:,:),         & ! init high-frequency stress                 (nproma,nblks_c) (m/s)^2
      &  phihf1(:,:),         & ! init high-frequency energy flux into ocean (nproma,nblks_c) (m/s)^2
      &  tauhf(:,:),          & ! high-frequency stress                      (nproma,nblks_c) (m/s)^2
      &  phihf(:,:),          & ! high-frequency energy flux into ocean      (nproma,nblks_c) (m/s)^2
      &  xlevtail(:,:),       & ! tail level                                 (nproma,nblks_c) (-)
      &  tauw(:,:),           & ! wave stress                                (nproma,nblks_c) (m/s)^2
      &  phiaw(:,:),          & ! energy flux from wind into waves integrated over the full frequency range  (nproma,nblks_c) (-)
      &  enh(:,:),            & ! nonlinear transfer function coefficients for shallow water 
      &  AF11(:),             & ! for discrete approximation of nonlinear transfer (nfreqs+4) (-)
      &  FKLAP(:), FKLAP1(:), & ! --//-- (nfreqs+4) (-)
      &  FKLAM(:), FKLAM1(:), & ! --//-- (nfreqs+4) (-)
      ! total waves
      &  emean(:,:),          & ! total energy                   (nproma,nblks_c) (m^2)
      &  emeanws(:,:),        & ! total wind sea input energy    (nproma,nblks_c) (m^2)
      &  femean(:,:),         & ! mean frequency energy          (nproma,nblks_c) (m^2)
      &  femeanws(:,:),       & ! windsea mean frequency energy  (nproma,nblks_c) (m^2)
      &  hs(:,:),             & ! total significant wave height  (nproma,nblks_c) (m)
      &  hs_dir(:,:),         & ! total mean wave direction      (nproma,nblks_c) (deg)
      &  tpp(:,:),            & ! total peak wave period         (nproma,nblks_c) (s)
      &  tmp(:,:),            & ! total mean wave period         (nproma,nblks_c) (s)
      &  tm1(:,:),            & ! total wave m1 period           (nproma,nblks_c) (s)
      &  tm2(:,:),            & ! total wave m2 period           (nproma,nblks_c) (s)
      &  ds(:,:),             & ! total directional wave spread  (nproma,nblks_c) (deg)
      &  hrms_frac(:,:),      & ! square ratio (Hrms / Hmax)**2  (nproma,nblks_c) (-)
      &  wbr_frac(:,:),       & ! fraction of breaking waves     (nproma,nblks_c) (-)
      ! wind sea
      &  emean_sea(:,:),      & ! wind sea energy                (nproma,nblks_c) (m^2)
      &  femean_sea(:,:),     & ! wind sea mean frequency energy (nproma,nblks_c) (m^2)
      &  f1mean_sea(:,:),     & ! wind sea mean frequency        (nproma,nblks_c) (Hz)
      &  hs_sea(:,:),         & ! sea significant wave height    (nproma,nblks_c) (m)
      &  hs_sea_dir(:,:),     & ! sea mean wave direction        (nproma,nblks_c) (deg)  
      &  pp_sea(:,:),         & ! sea peak period                (nproma,nblks_c) (s)  
      &  mp_sea(:,:),         & ! sea mean period                (nproma,nblks_c) (s)  
      &  m1_sea(:,:),         & ! sea m1-period                  (nproma,nblks_c) (s) 
      &  m2_sea(:,:),         & ! sea m2-period                  (nproma,nblks_c) (s)
      &  ds_sea(:,:),         & ! sea directional spreed         (nproma,nblks_c) (deg) 
      ! swell
      &  emean_swell(:,:),    & ! swell energy                   (nproma,nblks_c) (m^2)
      &  femean_swell(:,:),   & ! swell mean frequency energy    (nproma,nblks_c) (m^2)
      &  f1mean_swell(:,:),   & ! swell sea mean frequency       (nproma,nblks_c) (Hz)
      &  hs_swell(:,:),       & ! swell significant wave height  (nproma,nblks_c) (m)
      &  hs_swell_dir(:,:),   & ! swell mean wave direction      (nproma,nblks_c) (deg)
      &  pp_swell(:,:),       & ! swell peak period              (nproma,nblks_c) (s)
      &  mp_swell(:,:),       & ! swell mean period              (nproma,nblks_c) (s)
      &  m1_swell(:,:),       & ! swell m1-period                (nproma,nblks_c) (s)
      &  m2_swell(:,:),       & ! swell m2-period                (nproma,nblks_c) (s)
      &  ds_swell(:,:),       & ! swell directional spreed       (nproma,nblks_c) (deg)
      !
      &  drag(:,:),           & ! drag coefficient               (nproma,nblks_c) (-)
      &  tauwn(:,:),          & ! normalised wave stress         (nproma,nblks_c) (-)
      &  beta(:,:),           & ! Charnock parameter             (nproma,nblks_c) (-)
      &  mean_period(:,:),    & ! mean wave period = 1/femean    (nproma,nblks_c) (s)
      &  peak_period(:,:),    & ! peak wave period               (nproma,nblks_c) (s)
      &  u_stokes(:,:),       & ! U-component of surface Stokes drift (nproma,nblks_c) (m/s)
      &  v_stokes(:,:)        & ! V-component of surface Stokes drift (nproma,nblks_c) (m/s)
      &  => NULL()

    INTEGER, POINTER, CONTIGUOUS ::  &
      &  last_prog_freq_ind(:,:), & ! last frequency index of the prognostic range (nproma,nblks_c) (-)
      &  swell_mask(:,:),         & ! swell separation mask (nproma,nblks_c) (-)
      &  swell_mask_tr(:,:,:),    & ! swell separation mask (nproma,nblks_c,ntracer) (-)
      &  ikp(:), ikp1(:),         & ! for discrete approximation of nonlinear transfer (nfreqs+4) (-)
      &  ikm(:), ikm1(:),         & ! --//-- (nfreqs+4) (-)
      &  k1w(:,:), k2w(:,:),      & ! --//-- (ndirs, 2) (-)
      &  k11w(:,:), k21w(:,:),    & ! --//-- (ndirs, 2) (-) 
      &  ja1(:,:), ja2(:,:),      & ! --//-- (ndirs, 2) (-)
      &  non_lin_tr_ind(:,:,:,:)  & ! tracer index for nonlinear interaction (nfreqs+4,2,ndirs,8) (-)
      &  => NULL()

    TYPE(t_ptr_2d3d),   ALLOCATABLE :: gv_e_freq_ptr(:)  !< pointer array: one pointer for each frequence
    TYPE(t_ptr_2d3d),   ALLOCATABLE :: gv_c_freq_ptr(:)  !< pointer array: one pointer for each frequence
    TYPE(t_ptr_2d3d),   ALLOCATABLE :: wave_num_c_ptr(:) !< pointer array: one pointer for each frequence
    TYPE(t_ptr_2d3d),   ALLOCATABLE :: wave_num_e_ptr(:) !< pointer array: one pointer for each frequence
    TYPE(t_ptr_2d_int), ALLOCATABLE :: swmask_ptr(:)     !< pointer array: one pointer for each tracer
  END type t_wave_diag

  TYPE t_wave_state
    !array of prognostic states at different timelevels
    TYPE(t_wave_prog), ALLOCATABLE    :: prog(:)       !< shape: (timelevels)
    TYPE(t_wave_source)               :: source        !< source function state vector
    TYPE(t_wave_diag)                 :: diag
  END TYPE t_wave_state

  TYPE t_wave_state_lists
    ! array of prognostic state lists at different timelevels
    TYPE(t_var_list_ptr), ALLOCATABLE :: prog_list(:)  !< shape: (timelevels)
    TYPE(t_var_list_ptr)              :: source_list
    TYPE(t_var_list_ptr)              :: diag_list
  END TYPE t_wave_state_lists

END MODULE mo_wave_types
