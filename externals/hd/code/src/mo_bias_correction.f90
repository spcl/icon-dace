! mo_bias_correction.f90 - Module for bias correction of discharge at the river mouths
! 
! Copyright (C) 2023, Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
! SPDX-License-Identifier: Apache-2.0
! See ./LICENSES/ for license information
!
! Author: Stefan Hagemann
! Contact: <stefan.hagemann@hereon.de>
!_________________________________________

MODULE mo_bias_correction

  ! Vs. 1.0 - Nov. 2023 - Stefan Hagemann, Hereon

  USE mo_exception,      ONLY: message, message_text, finish
  USE mo_hydrology,      ONLY: grid_hd
  USE mo_io,             ONLY: io_open, io_close, io_read, io_write
  USE mo_kind,           ONLY: dp
  USE mo_mpi,            ONLY: p_io, p_parallel_io, p_parallel, p_bcast
  USE mo_netcdf,         ONLY: file_info, io_inq_varid, io_get_var_double
  USE mo_time_control,   ONLY: delta_time 

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: bc_init, bc_cleanup, bias_correction

!
! *** Bias Correction parameter
  REAL(dp), ALLOCATABLE :: perc_low(:,:)     ! Lower percentile iperc in model data
  REAL(dp), ALLOCATABLE :: perc_high(:,:)    ! Higher percentile 1-iperc in reference data
  REAL(dp), ALLOCATABLE :: factor_low(:,:)   ! Correction factor for model values below the iperc percentile 
  REAL(dp), ALLOCATABLE :: factor_high(:,:)  ! Correction factor for model values above the 1-iperc percentile 
  REAL(dp), ALLOCATABLE :: factor_mid(:,:)   ! Correction factor for model values inbetween the iperc and 1-iperc percentiles 
  REAL(dp), ALLOCATABLE :: factor_mean(:,:)  ! Correction factor for mean bias 
  REAL(dp), ALLOCATABLE, PUBLIC :: hd_bc_outflow(:,:) ! Bias corrected outflow
!
  REAL(dp), PARAMETER :: zeps = 1.E-10_dp
  REAL(dp), PARAMETER :: xmiss = -9999._dp

CONTAINS

! *******************************************************************************
  SUBROUTINE bc_init(dn_bcpara)
! *******************************************************************************
!
!   *** Initialize Bias Correction of Discharge within the HD model
!
    CHARACTER(LEN=*), INTENT(in) :: dn_bcpara    ! File with Bias correction parameters

    ! local variables
    INTEGER                :: read_status, inml, iunit

    ! Initialize BC memory
    CALL bc_init_memory(grid_hd%nlon, grid_hd%nlat)

    ! Read BC parameter from file
    CALL read_bc_parameter(dn_bcpara)

    WRITE (message_text,*) ' Bias correction initialized in bc_init'
    CALL message('bc_init', message_text)

  END SUBROUTINE bc_init

! *******************************************************************************
  SUBROUTINE bias_correction(ibc_type, forg)
! *******************************************************************************
!
!   *** Conduct bias correction for all river mouths and stores bias corrected outflow in hd_bc_outflow
!   *** Note that BC parameters are derived for daily discharges in m^3/s
!   *** Hence applying them to sub-daily discharges will not yield the same results as an 
!   *** offline application on daily discharges.
!     
!   *** del_smooth = Smoothing radius at percentiles borders
!
    INTEGER, INTENT(in) :: ibc_type   ! BC type: 0=None, 1 = Mean Bias, 2 = Low, Mid and High Biases            ! 
    REAL(dp), DIMENSION(:,:), INTENT(in)  :: forg   ! Uncorrected original discharges to be corrected
!
    REAL(dp) :: del_smooth = 0.05_dp               ! Smoothing radius at percentiles borders
    REAL(dp) :: xsm_upper 
    REAL(dp) :: xsm_lower
!
    xsm_upper = 1._dp + del_smooth 
    xsm_lower = 1._dp - del_smooth
!
!   *** Loop over river mouths 
    !#TODO  - take care on correct unit
!
!   ******* Conduct bias correction
    SELECT CASE(ibc_type) 
      CASE(1)            ! Mean Bias correction
        WHERE (ABS(factor_mean(:,:)-xmiss).GT.zeps)
          hd_bc_outflow = forg(:,:) * factor_mean(:,:) 
        ELSEWHERE
          hd_bc_outflow = forg(:,:)
        END WHERE
      CASE(2)            ! Bias correction for low, mid and high percentile ranges
        IF (del_smooth.LE.zeps) THEN
          WHERE (ABS(factor_mid(:,:)-xmiss).GT.zeps)
            WHERE (forg(:,:) .LE. perc_low(:,:))
              hd_bc_outflow = forg(:,:) * factor_low(:,:)
            ELSE WHERE (forg(:,:) .GE. perc_high(:,:))
              hd_bc_outflow = forg(:,:) * factor_high(:,:)
            ELSEWHERE
              hd_bc_outflow = forg(:,:) * factor_mid(:,:)
            END WHERE
          ELSEWHERE
            hd_bc_outflow = forg(:,:)
          END WHERE
        ELSE
          WHERE (ABS(factor_mid(:,:)-xmiss).GT.zeps)
            WHERE (forg(:,:) .LE. xsm_lower*perc_low(:,:))
              hd_bc_outflow = forg(:,:) * factor_low(:,:)
            ELSE WHERE (forg(:,:) .GE. xsm_upper*perc_high(:,:))
              hd_bc_outflow = forg(:,:) * factor_high(:,:)
            ! smooth transitions
            ELSE WHERE (forg(:,:) .LT. xsm_upper*perc_low(:,:))
              hd_bc_outflow = forg(:,:) * ( factor_low(:,:) + (factor_mid(:,:) - factor_low(:,:)) *  &
                     (forg(:,:) - xsm_lower*perc_low(:,:)) / &
                     (2._dp*del_smooth*perc_low(:,:)) )
            ELSE WHERE (forg(:,:) .GT. xsm_lower*perc_high(:,:))
              hd_bc_outflow = forg(:,:) * ( factor_mid(:,:) + (factor_high(:,:) - factor_mid(:,:)) *  &
                     (forg(:,:) - xsm_lower*perc_high(:,:)) / &
                     (2._dp*del_smooth*perc_high(:,:)) )
            ! mid range
            ELSEWHERE
              hd_bc_outflow = forg(:,:) * factor_mid(:,:)
            END WHERE
          ELSEWHERE
            hd_bc_outflow = forg(:,:)
          END WHERE
        ENDIF
      CASE DEFAULT
        hd_bc_outflow = forg(:,:)
    END SELECT

  END SUBROUTINE bias_correction
!
! *******************************************************************************
  SUBROUTINE bc_init_memory(nl, nb)
! *******************************************************************************
! 
!   *** Allocate BC Parameter arrays

    INTEGER,    INTENT(in)    :: nl      ! Number of HD longitudes
    INTEGER,    INTENT(in)    :: nb      ! Number of HD latitudes

    IF (p_parallel_io) THEN
      ALLOCATE (perc_low(nl,nb))          ; perc_low(:,:)    = 0.0_dp
      ALLOCATE (perc_high(nl,nb))         ; perc_high(:,:)   = 0.0_dp
      ALLOCATE (factor_low(nl,nb))        ; factor_low(:,:)  = 0.0_dp
      ALLOCATE (factor_mid(nl,nb))        ; factor_mid(:,:)  = 0.0_dp
      ALLOCATE (factor_high(nl,nb))       ; factor_high(:,:) = 0.0_dp
      ALLOCATE (factor_mean(nl,nb))       ; factor_mean(:,:) = 0.0_dp
      ALLOCATE (hd_bc_outflow(nl,nb))     ; hd_bc_outflow(:,:) = 0.0_dp
    ENDIF

  END SUBROUTINE bc_init_memory
!
! *******************************************************************************
  SUBROUTINE read_bc_parameter(dn_bcpara)
! *******************************************************************************
!
!   ******* Read Bias correction parameter file
!
    CHARACTER(LEN=*), INTENT(in) :: dn_bcpara
!
    TYPE (FILE_INFO)  :: parafile
    INTEGER :: fileID, nvarid, i
    INTEGER, PARAMETER  :: nvar = 6
    CHARACTER (LEN=20) :: cvar

    parafile%opened = .FALSE.
    CALL IO_open (dn_bcpara, parafile, IO_READ)
    CALL message('', '')
    WRITE (message_text,*) 'Reading Bias correction parameter from file ', TRIM(dn_bcpara)
    CALL message('read_bc_parameters', message_text)

    fileID = parafile%file_id

    DO i=1, nvar
      SELECT CASE(i)
        CASE(1) ; cvar = 'Factor_Mean' 
        CASE(2) ; cvar = 'Factor_Low'  
        CASE(3) ; cvar = 'Factor_Mid'  
        CASE(4) ; cvar = 'Factor_High' 
        CASE(5) ; cvar = 'Perc_Low'    
        CASE(6) ; cvar = 'Perc_High'   
      END SELECT 
      CALL IO_inq_varid (fileID, cvar, nvarid)
      SELECT CASE(i)
        CASE(1) ; CALL IO_get_var_double (fileID, nvarid, factor_mean)
        CASE(2) ; CALL IO_get_var_double (fileID, nvarid, factor_low)
        CASE(3) ; CALL IO_get_var_double (fileID, nvarid, factor_mid)
        CASE(4) ; CALL IO_get_var_double (fileID, nvarid, factor_high)
        CASE(5) ; CALL IO_get_var_double (fileID, nvarid, perc_low)  
        CASE(6) ; CALL IO_get_var_double (fileID, nvarid, perc_high)
      END SELECT 
    ENDDO

    CALL IO_close(parafile)
!
  END SUBROUTINE read_bc_parameter
!
! *******************************************************************************
  SUBROUTINE bc_cleanup
! *******************************************************************************
!
!   Free memory used for bias correction

    IF (p_parallel_io) THEN
      DEALLOCATE  (perc_low)
      DEALLOCATE  (perc_high)
      DEALLOCATE  (factor_low)
      DEALLOCATE  (factor_mid)
      DEALLOCATE  (factor_high)
      DEALLOCATE  (factor_mean)
      DEALLOCATE  (hd_bc_outflow)
    ENDIF

  END SUBROUTINE bc_cleanup
!
END MODULE mo_bias_correction


