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

! Main program for the ICON atmospheric model

!-----------------------------
#include "omp_definitions.inc"
!-----------------------------

MODULE mo_test_jitter

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, message_text, finish
  USE mo_mpi,                 ONLY: work_mpi_barrier, my_process_is_stdio
  USE mo_timer,               ONLY: new_timer, timer_start, timer_stop, &
    & print_timer, cleanup_timer

  USE mo_master_control,      ONLY: get_my_process_name, get_my_model_no
  USE mo_icon_testbed_config, ONLY: testbed_iterations, calculate_iterations, &
    & no_of_blocks, no_of_layers

  USE mo_parallel_nml,       ONLY: read_parallel_namelist
  USE mo_parallel_config,    ONLY: nproma

!-------------------------------------------------------------------------

IMPLICIT NONE

PRIVATE

PUBLIC :: test_jitter

CONTAINS
  
  !-------------------------------------------------------------------------
  !>
  !!
  SUBROUTINE test_jitter(namelist_filename,shr_namelist_filename)

    CHARACTER(LEN=*), INTENT(in) :: namelist_filename
    CHARACTER(LEN=*), INTENT(in) :: shr_namelist_filename
        
    INTEGER ::  timer_barrier_init, iter
        
    CHARACTER(*), PARAMETER :: method_name = "mo_test_jitter:test_jitter"


    !---------------------------------------------------------------------
    CALL read_parallel_namelist(namelist_filename)
    
    CALL message(" ---------------- ", method_name)
    WRITE(message_text,*) "testbed_iterations=", testbed_iterations
    CALL message(" -- ", message_text)
    WRITE(message_text,*) "calculate_iterations=", calculate_iterations
    CALL message(" -- ", message_text)
    WRITE(message_text,*) "nproma=", nproma, " layers=", no_of_layers, " blocks=", no_of_blocks
    CALL message(" -- ", message_text)
   !-------------------------------------------------------------------------
    timer_barrier_init  = new_timer("mpi_barrier_init")
    CALL work_mpi_barrier()
    
    !---------------------------------------------------------------------
    ! call some barriers to see how much time it takes
    DO iter=1, testbed_iterations
      CALL timer_start(timer_barrier_init)
      CALL work_mpi_barrier()    
      CALL timer_stop(timer_barrier_init)
    ENDDO
    !---------------------------------------------------------------------
    CALL test_jitter_iter()
    !---------------------------------------------------------------------
  
  END SUBROUTINE test_jitter
  !---------------------------------------------------------------------
    
  !---------------------------------------------------------------------
  !>
  SUBROUTINE test_jitter_iter()

    ! 3D variables
    REAL(wp), DIMENSION(nproma,no_of_layers,no_of_blocks) :: a, b, c
    REAL(wp) :: suma, sumb, sumc
    
    INTEGER :: timer_calculate, timer_barrier
    
    INTEGER :: i, j, k, iter, calculate
    
    CHARACTER(*), PARAMETER :: method_name = "mo_test_jitter:test_jitter_iter"
    
    timer_barrier  = new_timer("mpi_barrier")
    timer_calculate  = new_timer("calculate")
!$OMP PARALLEL
!$OMP DO PRIVATE(i,k,j) ICON_OMP_DEFAULT_SCHEDULE
    DO i = 1, no_of_blocks
      DO k = 1, no_of_layers
        DO j = 1, nproma
          a(j,k,i) = 0.0_wp
          b(j,k,i) = real(k,wp)
          c(j,k,i) = real(k+i,wp)
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
    !---------------------------------------------------------------------
    CALL work_mpi_barrier()    
    !---------------------------------------------------------------------

    DO iter=1, testbed_iterations
      !---------------------------------------------------------------------
      ! do some calculations
      CALL timer_start(timer_calculate)
!$OMP PARALLEL
      DO calculate=1,calculate_iterations

!$OMP DO PRIVATE(i,k,j) ICON_OMP_DEFAULT_SCHEDULE
        DO i = 1, no_of_blocks
          DO k = 1, no_of_layers
            DO j = 1, nproma
              a(j,k,i) = b(j,k,i) / c(j,k,i)
            ENDDO
          ENDDO
        ENDDO
!$OMP END DO

!$OMP DO PRIVATE(i,k,j) ICON_OMP_DEFAULT_SCHEDULE
        DO i = 1, no_of_blocks
          DO k = 1, no_of_layers
            DO j = 1, nproma
              b(j,k,i) = c(j,k,i) - a(j,k,i)
            ENDDO
          ENDDO
        ENDDO
!$OMP END DO
                
!$OMP DO PRIVATE(i,k,j) ICON_OMP_DEFAULT_SCHEDULE
        DO i = 1, no_of_blocks
          DO k = 1, no_of_layers
            DO j = 1, nproma
              c(j,k,i) = MAX(a(j,k,i) * b(j,k,i), ABS(b(j,k,i)-a(j,k,i)))
            ENDDO
          ENDDO
        ENDDO
!$OMP END DO


      ENDDO !calculate=1,calculate_iterations
!$OMP END PARALLEL
            
      CALL timer_stop(timer_calculate)
!       write(0,*) c(nproma,no_of_layers,no_of_blocks)
      !---------------------------------------------------------------------
      CALL timer_start(timer_barrier)
      CALL work_mpi_barrier()    
      CALL timer_stop(timer_barrier)
      CALL work_mpi_barrier()    
      !---------------------------------------------------------------------
          
    ENDDO !iter=1, testbed_iterations
             
    !---------------------------------------------------------------------
    ! print something to avoid optimization misfortunes
    suma = SUM(a(:,:,:))
    sumb = SUM(b(:,:,:))
    sumc = SUM(c(:,:,:))
    write(0,*) "sums=", suma, sumb, sumc
    !---------------------------------------------------------------------
    
    !---------------------------------------------------------------------
    ! print the timers
    CALL print_timer()
    CALL cleanup_timer(timer_calculate)
    CALL cleanup_timer(timer_barrier)
    !---------------------------------------------------------------------

  END SUBROUTINE test_jitter_iter
  !-------------------------------------------------------------------------


END MODULE mo_test_jitter

