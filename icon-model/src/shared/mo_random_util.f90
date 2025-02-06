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

! This module contains the subroutine that Adds random perturbation
! to the normal wind for the Non-Hyd core

MODULE mo_random_util

  USE mo_kind,                     ONLY: wp, i4, i8
  USE mo_impl_constants,           ONLY: SUCCESS, ON_CELLS, ON_EDGES, ON_VERTICES
  USE mo_exception,                ONLY: message, finish, message_text
  USE mo_grid_subset,              ONLY: t_subset_range, get_index_range
  USE mo_math_types,               ONLY: t_geographical_coordinates
  USE mo_random_number_generators, ONLY: initialize_random_number_generator, generate_uniform_random_number
 ! USE mo_mpi,                 ONLY: get_my_global_mpi_id, global_mpi_barrier

  IMPLICIT NONE
  PRIVATE
 
  PUBLIC :: add_random_noise_global, add_random_noise_3d, add_random_noise_2d

CONTAINS

  !-----------------------------------------------------------------------
  !>
  !! Add random perturbation using global index
  !! for seeds, thus guaranteeing reproducability
  !!
  SUBROUTINE add_random_noise_global( &
    & in_subset,                 & ! in
    & in_var,                    & ! inout
    & start_level, end_level,    & ! in
    & noise_scale,               & ! in
    & global_vertical_seed,      & ! input, OPTIONAL
    & add_seed,                  & ! input, OPTIONAL
    & debug)                       ! input, OPTIONAL

    TYPE(t_subset_range)    :: in_subset
    REAL(wp), INTENT(INOUT) :: in_var(:,:,:) ! variable to perturb
    INTEGER, INTENT(IN)     :: start_level, end_level
    REAL(wp), INTENT(IN)    :: noise_scale ! magnitude of the perturbation
    INTEGER, INTENT(IN), TARGET, OPTIONAL  :: global_vertical_seed(:)
    INTEGER, INTENT(IN), OPTIONAL    :: add_seed ! add this number to the fisrt seed element to obtain the rest of them
    LOGICAL, INTENT(IN), OPTIONAL :: debug


    ! LOCAL VARIABLES
    ! INTEGER :: DateTimeArray(8)    ! Holds the date and time
    INTEGER :: status
    INTEGER :: level, block, idx, start_idx, end_idx

    INTEGER :: seed_size, js, seed_trigger
    INTEGER, ALLOCATABLE :: seed_array(:)

    REAL(wp), ALLOCATABLE :: noise_1D(:), noise_3D(:,:,:)
    INTEGER, POINTER  :: vertical_seed(:)
    INTEGER    :: add_thisSeed
    LOGICAL :: mydebug

    CHARACTER(len=*), PARAMETER ::  &
      &  method_name = 'mo_random_util:add_random_noise_global'

    !-----
    mydebug = .FALSE.
    IF (PRESENT(debug)) mydebug = debug

    IF (mydebug) THEN
      CALL message(TRIM(method_name),'=========== generating random noise based on global index =============')
    END IF

    !-----------------------------------------------------------
    ! 1. prepare memory for the seed
    !-----------------------------------------------------------

    CALL RANDOM_SEED(SIZE=seed_size)
    WRITE(message_text,*) 'The size of the intrinsic seed is', seed_size
    IF (mydebug) THEN
      CALL message(method_name,TRIM(message_text))
    END IF

    ALLOCATE( seed_array(seed_size), STAT=status)
    IF(status/=SUCCESS)THEN
      CALL finish(method_name,'allocation of seed_array failed')
    ENDIF
 
    IF (PRESENT(global_vertical_seed)) THEN
      vertical_seed => global_vertical_seed
    ELSE
      SELECT CASE (in_subset%entity_location)
      CASE (ON_CELLS )
        vertical_seed => in_subset%patch%cells%decomp_info%glb_index
      CASE (ON_EDGES )
        vertical_seed => in_subset%patch%edges%decomp_info%glb_index
      CASE (ON_VERTICES )
        vertical_seed => in_subset%patch%verts%decomp_info%glb_index
      CASE default
        CALL finish(method_name, "Unknown in_subset%entity_location")
      END SELECT
    END IF
    add_thisSeed = 7
    IF (PRESENT(add_seed)) add_thisSeed = add_seed


    ! use the global index as a seed in order to ensure p_test works
    ! this implies a random sequence for each column
    ALLOCATE(noise_1D(start_level:end_level), STAT=status)

    DO block = in_subset%start_block, in_subset%end_block
      CALL get_index_range( in_subset, block, start_idx, end_idx)
      DO idx = start_idx, end_idx
        seed_trigger = ABS(vertical_seed( (block-1) * in_subset%block_size + idx) )
        DO js = 0, seed_size-1
           seed_array(js+1) = seed_trigger + js * add_thisSeed
        ENDDO
        ! Since this routine is not called from OpenMP, no need to change the seed
        ! WARNING: changing the seed produced LOW-QUALITY and CORRELATED sequences
        ! WS 2019-11-15: reenabled in order to pass Mistral atm_rce_les tests
        !                but Dmitry's comment above has to be addressed...
        CALL RANDOM_SEED( PUT=seed_array )
        CALL RANDOM_NUMBER( noise_1D(start_level:end_level))

!          WRITE(0,*) get_my_global_mpi_id(), ":", vertical_seed( (block-1) * in_subset%block_size + idx), &
!            " seeds:", seed_array(:), " noise=", noise_1D(start_level:end_level)

        DO level = start_level, end_level
          in_var(idx,level,block) = in_var(idx,level,block) + (noise_1D(level) - 0.5_wp) * noise_scale
        ENDDO

      ENDDO ! idx = start_idx, end_idx

    ENDDO ! block = in_subset%start_block, in_subset%end_block

    DEALLOCATE( noise_1D )

 !   CALL global_mpi_barrier()
 !   CALL finish(method_name, " global_vertical_seed")

  END SUBROUTINE add_random_noise_global

  SUBROUTINE add_random_noise_3d(subset, coordinates, amplitude, field, seed_in)

    TYPE(t_subset_range), INTENT(IN) :: subset
    TYPE(t_geographical_coordinates), ALLOCATABLE, DIMENSION(:,:), INTENT(IN) :: coordinates

    REAL(wp), INTENT(IN) :: amplitude
    INTEGER(i8), INTENT(IN) :: seed_in
    REAL(wp), INTENT(INOUT) :: field(:,:,:)

    INTEGER :: rng_state, iseed, jseed, start_idx, end_idx, jl, jk, jb

    iseed = IEOR(INT(seed_in/2_i8**32_i8, i4), INT(seed_in, i4))

    DO jb=subset%start_block,subset%end_block
      CALL get_index_range( subset, jb, start_idx, end_idx)
      DO jk=1,SIZE(field,2)
        DO jl=start_idx,end_idx
          jseed = IEOR(TRANSFER(coordinates(jl,jb)%lat,jseed), TRANSFER(coordinates(jl,jb)%lon,jseed)*2**9)
          jseed = IEOR(jseed, jk*2**7)
          rng_state = initialize_random_number_generator(iseed, jseed)
          field(jl,jk,jb) = field(jl,jk,jb) * ((generate_uniform_random_number(rng_state) * 2._wp - 1._wp) * amplitude + 1._wp)
        ENDDO !jl
      ENDDO !jk
    ENDDO !jb

  END SUBROUTINE add_random_noise_3d

  SUBROUTINE add_random_noise_2d(subset, coordinates, amplitude, field, seed_in)

    TYPE(t_subset_range), INTENT(IN) :: subset
    TYPE(t_geographical_coordinates), ALLOCATABLE, DIMENSION(:,:), INTENT(IN) :: coordinates

    REAL(wp), INTENT(IN) :: amplitude
    INTEGER(i8), INTENT(IN) :: seed_in
    REAL(wp), INTENT(INOUT) :: field(:,:)

    INTEGER :: rng_state, iseed, jseed, start_idx, end_idx, jl, jb

    iseed = IEOR(INT(seed_in/2_i8**32_i8, i4), INT(seed_in, i4))

    DO jb=subset%start_block,subset%end_block
      CALL get_index_range( subset, jb, start_idx, end_idx)
        DO jl=start_idx,end_idx
          jseed = IEOR(TRANSFER(coordinates(jl,jb)%lat,jseed), TRANSFER(coordinates(jl,jb)%lon,jseed)*2**9)
          rng_state = initialize_random_number_generator(iseed, jseed)
          field(jl,jb) = field(jl,jb) * ((generate_uniform_random_number(rng_state) * 2._wp - 1._wp) * amplitude + 1._wp)
        ENDDO !jl
    ENDDO !jb

  END SUBROUTINE add_random_noise_2d

END MODULE mo_random_util
