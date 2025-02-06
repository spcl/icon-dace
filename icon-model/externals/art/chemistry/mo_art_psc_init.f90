!
! mo_art_psc_init
! This module provides routines for initialisation of PSCs
! It has to be separated from mo_art_psc_state because then creation of a
! standalone is much easier
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

MODULE mo_art_psc_init

! ICON
  USE mo_kind,                          ONLY: wp
  USE mo_exception,                     ONLY: finish 
  USE mo_var_list,                      ONLY: t_var_list_ptr
  USE mo_var,                           ONLY: t_var
  USE mo_tracer_metadata_types,         ONLY: t_chem_meta
  USE mo_var_metadata_types,            ONLY: t_var_metadata_dynamic, t_var_metadata
  USE mo_art_config,                    ONLY: art_config
  USE mo_impl_constants,                ONLY: SUCCESS
! ART
  USE mo_art_data,                      ONLY: p_art_data
  USE mo_art_chem_data,                 ONLY: t_art_chem_indices
  USE mo_art_atmo_data,                 ONLY: t_art_atmo
  USE mo_art_psc_types,                 ONLY: t_art_psc
  USE messy_mecca_kpp_global,           ONLY: IHS_MAX
  
  IMPLICIT NONE
  
  PRIVATE

  PUBLIC :: art_psc_init_arrays
  PUBLIC :: art_psc_set_number_bins

  
CONTAINS
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE art_psc_set_number_bins(PSC,jg)
!<
! SUBROUTINE art_psc_init_arrays
! Get the number of size bins given in the  chemtracer XML
! Part of Module: mo_art_psc_state
! Author: Michael Weimer, KIT
! Initial Release: 2020-06-19
! Modifications: 
!>
  IMPLICIT NONE
  TYPE(t_art_psc), INTENT(inout) :: &
    &  PSC                         !< PSC meta data structure
  INTEGER, INTENT(in) :: &
    &  jg                          !< patch id
  ! local variables
  CHARACTER(LEN = 2) :: NSB_str    !< NSB converted to character
  INTEGER ::   &
    &  ierror, &                   !< index if tracer of nth size bin is found
    &  iTR                         !< index of nth size bin (if found)
    
  PSC%NSB = 0
  ierror = SUCCESS

  DO WHILE(ierror == SUCCESS)
    PSC%NSB = PSC%NSB + 1
    WRITE(NSB_str,'(I2.2)') PSC%NSB

    CALL p_art_data(jg)%dict_tracer%get('TRNAT_bin'//NSB_str,iTR,ierror)
  END DO
  
  IF (PSC%NSB > 1) THEN
    PSC%NSB = PSC%NSB - 1
  ELSE
    PSC%NSB = 1
  END IF

END SUBROUTINE art_psc_set_number_bins
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE art_psc_init_arrays(PSC,jg,tracer,p_prog_list)
!<
! SUBROUTINE art_psc_init_arrays
! This subroutine allocates the arrays for the psc structure
! Part of Module: mo_art_psc_state
! Author: Michael Weimer, KIT
! Initial Release: 2016-11-21
! Modifications: 
! 2017-06-09: Michael Weimer, KIT
! - added initialisation of kinetic NAT parametrisation structure
!>
  IMPLICIT NONE
  TYPE(t_art_psc), INTENT(inout) :: &
    &  PSC                         !< PSC meta data structure
  INTEGER, INTENT(in) :: &
    &  jg                          !< patch id
  REAL(wp), POINTER, INTENT(inout),OPTIONAL :: &
    &  tracer(:,:,:,:)             !< Tracer mixing ratios [kg kg-1]
  TYPE(t_var_list_ptr),INTENT(in) :: &
    &  p_prog_list                 !< current list: prognostic
  ! local variables
  INTEGER ::   &
    &  ierror, &                   !< index if tracer of nth size bin is found
    &  iTR,    &                   !< index of nth size bin (if found)
    &  n,iv                        !< loop index over bins and list-elements
  INTEGER, POINTER :: &
    &  jsp                         !< tracer index of current tracer in p_prog_list

  TYPE(t_art_atmo), POINTER :: &
    &  art_atmo

  TYPE(t_art_chem_indices), POINTER :: &
    &  art_indices

  TYPE(t_var_metadata_dynamic), POINTER :: info_dyn  !< returns reference to tracer
                                                     !  metadata of current element
  TYPE(t_var_metadata), POINTER :: info              !< returns reference to tracer
                                                     !  metadata of current element

  CHARACTER(LEN = 2) :: NSB_str, n_str  !< NSB and n converted to character
  CHARACTER(LEN = 11) ::   &
    & total_max_dens_str, sum_densities_str !< limits of max number density in
                                            !  bin and sum of the max no
                                            !  densities as character 
                                            !  (for error checking)

  art_atmo  => p_art_data(jg)%atmo
  art_indices => p_art_data(jg)%chem%indices

  ! get the tracer indices if they are present
  IF (art_config(jg)%lart_chemtracer) THEN
    PSC%iTRHNO3   = art_indices%iTRHNO3
    PSC%iTRH2SO4  = art_indices%iTRH2SO4
    PSC%iTRN2O5   = art_indices%iTRN2O5
    PSC%iTRHCl    = art_indices%iTRHCl
    PSC%iTRHBr    = art_indices%iTRHBr
    PSC%iTRHOCl   = art_indices%iTRHOCl
    PSC%iTRHOBr   = art_indices%iTRHOBr
    PSC%iTRClONO2 = art_indices%iTRClONO2
    PSC%iTRBrONO2 = art_indices%iTRBrONO2
  END IF

  IF (art_config(jg)%lart_mecca) THEN
    CALL p_art_data(jg)%dict_tracer%get('HNO3',PSC%iTRHNO3,ierror)
      IF (ierror /= SUCCESS) PSC%iTRHNO3 = 0
    CALL p_art_data(jg)%dict_tracer%get('H2SO4',PSC%iTRH2SO4,ierror)
      IF (ierror /= SUCCESS) PSC%iTRH2SO4 = 0

    ! As soon as one of HNO3 or H2SO4 is not available in the MECCA chemistry mechanism,
    ! use the chemtracers. Below it is checked if both of them are there
    IF (ANY((/ PSC%iTRHNO3, PSC%iTRH2SO4 /) == 0)) THEN
      PSC%iTRHNO3   = art_indices%iTRHNO3
      PSC%iTRH2SO4  = art_indices%iTRH2SO4
      PSC%iTRN2O5   = art_indices%iTRN2O5
      PSC%iTRHCl    = art_indices%iTRHCl
      PSC%iTRHBr    = art_indices%iTRHBr
      PSC%iTRHOCl   = art_indices%iTRHOCl
      PSC%iTRHOBr   = art_indices%iTRHOBr
      PSC%iTRClONO2 = art_indices%iTRClONO2
      PSC%iTRBrONO2 = art_indices%iTRBrONO2
    ELSE
      CALL p_art_data(jg)%dict_tracer%get('N2O5',PSC%iTRN2O5,ierror)
        IF (ierror /= SUCCESS) PSC%iTRN2O5 = 0
      CALL p_art_data(jg)%dict_tracer%get('HCl',PSC%iTRHCl,ierror)
        IF (ierror /= SUCCESS) PSC%iTRHCl = 0
      CALL p_art_data(jg)%dict_tracer%get('HBr',PSC%iTRHBr,ierror)
        IF (ierror /= SUCCESS) PSC%iTRHBr = 0
      CALL p_art_data(jg)%dict_tracer%get('HOCl',PSC%iTRHOCl,ierror)
        IF (ierror /= SUCCESS) PSC%iTRHOCl = 0
      CALL p_art_data(jg)%dict_tracer%get('HOBr',PSC%iTRHOBr,ierror)
        IF (ierror /= SUCCESS) PSC%iTRHOBr = 0
      CALL p_art_data(jg)%dict_tracer%get('ClNO3',PSC%iTRClONO2,ierror)
        IF (ierror /= SUCCESS) PSC%iTRClONO2 = 0
      CALL p_art_data(jg)%dict_tracer%get('BrNO3',PSC%iTRBrONO2,ierror)
        IF (ierror /= SUCCESS) PSC%iTRBrONO2 = 0
    END IF
  END IF


  IF (PSC%iTRH2SO4 == 0) THEN
    CALL finish('mo_art_psc_state:art_psc_init_arrays',  &
           &    'Could not find tracer TRH2SO4 or H2SO4, although lart_psc is set to true.')
  END IF

  IF (PSC%iTRHNO3 == 0) THEN
    CALL finish('mo_art_psc_state:art_psc_init_arrays',  &
            &   'Could not find tracer TRHNO3 or HNO3, although lart_psc is set to true.')
  END IF

  ! upper and lowermost level where PSCs are computed
  PSC%kstart = 1
  PSC%kend = art_atmo%nlev


  IF (.NOT. PRESENT(tracer)) THEN
    CALL finish('mo_art_psc_state:art_psc_init_arrays',  &
      &         'Parameter tracer not present.')
  END IF

  ALLOCATE(PSC%k_het(art_atmo%nproma,art_atmo%nlev,art_atmo%nblks,IHS_MAX))

  ALLOCATE(PSC%HNO3_Nconc_g(art_atmo%nproma,art_atmo%nlev))
  ! PSC chemistry arrays
  ALLOCATE(PSC%v_th(art_atmo%nproma,art_atmo%nlev))
  ALLOCATE(PSC%mean_free_path(art_atmo%nproma,art_atmo%nlev))
  ALLOCATE(PSC%Nconc_higher(art_atmo%nproma,art_atmo%nlev))
  ! STS arrays
  ALLOCATE(PSC%iexception(art_atmo%nproma,art_atmo%nlev))
  ALLOCATE(PSC%sqrtt(art_atmo%nproma,art_atmo%nlev))
  ALLOCATE(PSC%pw(art_atmo%nproma,art_atmo%nlev))
  ALLOCATE(PSC%logpw(art_atmo%nproma,art_atmo%nlev))
  ALLOCATE(PSC%ns(art_atmo%nproma,art_atmo%nlev))
  ALLOCATE(PSC%pn0(art_atmo%nproma,art_atmo%nlev))
  ALLOCATE(PSC%phocl0(art_atmo%nproma,art_atmo%nlev))
  ALLOCATE(PSC%phobr0(art_atmo%nproma,art_atmo%nlev))
  ALLOCATE(PSC%phcl0(art_atmo%nproma,art_atmo%nlev))
  ALLOCATE(PSC%phbr0(art_atmo%nproma,art_atmo%nlev))
  ALLOCATE(PSC%parthno3(art_atmo%nproma,art_atmo%nlev))
  ALLOCATE(PSC%parthcl(art_atmo%nproma,art_atmo%nlev))
  ALLOCATE(PSC%parthbr(art_atmo%nproma,art_atmo%nlev))
  ALLOCATE(PSC%hhocl(art_atmo%nproma,art_atmo%nlev))
  ALLOCATE(PSC%hhobr(art_atmo%nproma,art_atmo%nlev))
  ALLOCATE(PSC%hhbr(art_atmo%nproma,art_atmo%nlev))
  ALLOCATE(PSC%mn(art_atmo%nproma,art_atmo%nlev))
  ALLOCATE(PSC%ms(art_atmo%nproma,art_atmo%nlev))
  ALLOCATE(PSC%wn(art_atmo%nproma,art_atmo%nlev))
  ALLOCATE(PSC%ws(art_atmo%nproma,art_atmo%nlev))
  ALLOCATE(PSC%whocl(art_atmo%nproma,art_atmo%nlev))
  ALLOCATE(PSC%whobr(art_atmo%nproma,art_atmo%nlev))
  ALLOCATE(PSC%wcl(art_atmo%nproma,art_atmo%nlev))
  ALLOCATE(PSC%tice(art_atmo%nproma,art_atmo%nlev))
  ALLOCATE(PSC%density(art_atmo%nproma,art_atmo%nlev))
  ALLOCATE(PSC%diff(art_atmo%nproma,art_atmo%nlev))
  ALLOCATE(PSC%wnen(art_atmo%nproma,art_atmo%nlev))
  ALLOCATE(PSC%internal_temp(art_atmo%nproma,art_atmo%nlev))
  

  ! initialise arrays with zero
  PSC%liqsur = 0._wp
  PSC%cgaml = 0._wp
  PSC%k_het = 0._wp
  PSC%HNO3_Nconc_g = 0._wp
  PSC%HNO3_Nconc_l = 0._wp
  PSC%HNO3_Nconc_s = 0._wp
  PSC%ice_vmr_Marti = 0._wp
  PSC%v_th = 0._wp
  PSC%mean_free_path = 0._wp
  PSC%Nconc_higher = 0._wp
  PSC%iexception = 0
  PSC%sqrtt = 0._wp
  PSC%pw = 0._wp
  PSC%logpw = 0._wp
  PSC%ns = 0._wp
  PSC%pn0 = 0._wp
  PSC%phocl0 = 0._wp
  PSC%phobr0 = 0._wp
  PSC%phcl0 = 0._wp
  PSC%phbr0 = 0._wp
  PSC%parthno3 = 0._wp
  PSC%parthcl = 0._wp
  PSC%parthbr = 0._wp
  PSC%hhocl = 0._wp
  PSC%hhobr = 0._wp
  PSC%hhbr = 0._wp
  PSC%mn = 0._wp
  PSC%ms = 0._wp
  PSC%wn = 0._wp
  PSC%ws = 0._wp
  PSC%whocl = 0._wp
  PSC%whobr = 0._wp
  PSC%wcl = 0._wp
  PSC%tice = 0._wp
  PSC%density = 0._wp
  PSC%diff = 0._wp
  PSC%wnen = 0._wp
  PSC%internal_temp = 0._wp


  ! if number of size bins is greater than 1, kinetic NAT parametrisation is
  ! used else the diagnostic parametrisation
  IF (PSC%NSB > 1) THEN
    ALLOCATE(PSC%rbin_min(PSC%NSB))
    ALLOCATE(PSC%rbin_max(PSC%NSB))
    ALLOCATE(PSC%rbin_av(PSC%NSB))
    ALLOCATE(PSC%no_density_limit(PSC%NSB))
    ALLOCATE(PSC%tracer_indices_bins(PSC%NSB))

    ! get the tracer indices of NAT bins
    DO n=1,PSC%NSB
      WRITE(NSB_str,'(I2.2)') n
      CALL p_art_data(jg)%dict_tracer%get('TRNAT_bin'//NSB_str,iTR,ierror)
      PSC%tracer_indices_bins(n) = iTR
    !  tracer(:,:,:,iTR) = 0.0_wp  ! This here is incompatible with restarts.
                                   ! The initialization (if needed) needs to 
                                   ! be done somewhere else
    END DO

    ! ----------------------------------
    ! --- read meta information minimum radius and maximum number density of
    ! --- each NAT bin
    ! ----------------------------------
    ! --- Start DO-loop over elements in list
    ! ----------------------------------

    DO iv = 1, p_prog_list%p%nvars

      ! ----------------------------------
      ! --- Get meta data of current element and assure that current element is tracer
      ! ----------------------------------

      info_dyn=>p_prog_list%p%vl(iv)%p%info_dyn
      info=>p_prog_list%p%vl(iv)%p%info

      IF (info_dyn%tracer%lis_tracer) THEN

        ! ----------------------------------
        ! --- Retrieve running index
        ! ----------------------------------

        jsp=>info%ncontained

        IF (ANY(jsp == PSC%tracer_indices_bins(:))) THEN
          SELECT TYPE(meta => info_dyn%tracer)
            CLASS IS (t_chem_meta)
              n = 1
              DO WHILE ((n <= PSC%NSB) .AND. (jsp /= PSC%tracer_indices_bins(n)))
                n = n + 1
              END DO

              IF (n > PSC%NSB) THEN
                CALL finish('mo_art_psc_state:art_psc_init_arrays', & 
                       &    'could not find tracer '//TRIM(meta%name))
              END IF
 
              CALL meta%opt_meta%get('rbin_min',PSC%rbin_min(n),ierror)
              IF (ierror /= 0) THEN
                CALL finish('mo_art_psc_state:art_psc_init_arrays', &
                        &   'could not find meta data rbin_min for tracer' &
                        & //TRIM(meta%name))
              END IF
              CALL meta%opt_meta%get('no_density_limit',PSC%no_density_limit(n),ierror)
              IF (ierror /= 0) THEN
                CALL finish('mo_art_psc_state:art_psc_init_arrays',                &
                        &   'could not find meta data no_density_limit for tracer' &
                        & //TRIM(meta%name))
              END IF
          END SELECT
        END IF
      END IF
   
    END DO
    
    ! calculate maximum and average bin size
    DO n = 1,PSC%NSB
      IF (n < PSC%NSB) THEN
        IF (PSC%rbin_min(n+1) <= PSC%rbin_min(n)) THEN
          WRITE(n_str,'(I2)') n
          CALL finish('mo_art_psc_state:art_psc_init_arrays', &
                  &   'maximum is lower equal than minimum radius of bin with number '//n_str)
        END IF

        PSC%rbin_max(n) = PSC%rbin_min(n+1)
      ELSE
        PSC%rbin_max(n) = 25.e-6_wp ! should be 25 mu m
      END IF

      PSC%rbin_av(n) = (PSC%rbin_max(n) + PSC%rbin_min(n)) / 2._wp
    END DO

    ! total number density must be equal to total_max_density
    IF ((SUM(PSC%no_density_limit(1:PSC%NSB-1)) > PSC%total_max_density*1.0001) &
      .OR. (SUM(PSC%no_density_limit(1:PSC%NSB-1)) < PSC%total_max_density*0.9999)) THEN
  
      WRITE(sum_densities_str, '(ES11.4)') SUM(PSC%no_density_limit(1:PSC%NSB-1))
      WRITE(total_max_dens_str,'(ES11.4)') PSC%total_max_density
      CALL finish('mo_art_psc_state:art_psc_init_arrays',                           &
             &    'The sum of all no_density_limits must be equal to the given '    &
             &  //'total number density of '//total_max_dens_str//'. Current sum: ' &
             &  // sum_densities_str)
    END IF
    ALLOCATE(PSC%v_sed_NAT(art_atmo%nproma,PSC%kstart-1:PSC%kend+1,art_atmo%nblks,PSC%NSB))

  ELSE
    ! in case of the diagnostic NAT parametrisation just allocate the diagnostic
    ! arrays, use special condition: number of size bins = 1
    PSC%NSB = 1

    ALLOCATE(PSC%v_sed_NAT(art_atmo%nproma,PSC%kstart-1:PSC%kend+1,art_atmo%nblks,PSC%NSB))

  END IF !NSB > 1

  ! initialise arrays with zero
  PSC%dens_NAT = 0._wp
  PSC%radius_NAT = 0._wp
  PSC%NAT_sedi_rel_diff = 0._wp


END SUBROUTINE art_psc_init_arrays
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
END MODULE mo_art_psc_init
