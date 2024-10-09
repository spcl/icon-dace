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

! @brief configuration setup for z/i/p-level output
!
! configuration setup for z/i/p-level output

MODULE mo_nh_pzlev_config

  USE mo_kind,               ONLY: wp
  USE mo_impl_constants,     ONLY: SUCCESS, MAX_CHAR_LENGTH, max_dom
  USE mo_math_utilities,     ONLY: t_value_set, deallocate_set
  USE mo_exception,          ONLY: finish
  USE mo_util_sort,          ONLY: quicksort
  USE mo_mpi,                ONLY: my_process_is_stdio

  IMPLICIT NONE
  PUBLIC



  !!--------------------------------------------------------------------------
  !! Basic configuration setup for z/p-level output
  !!--------------------------------------------------------------------------
  TYPE :: t_nh_pzlev_config

    ! namelist variables
    !
    TYPE (t_value_set) :: zlevels    !< zlevel heights [m] 
    TYPE (t_value_set) :: plevels    !< plevel heights [Pa] 
    TYPE (t_value_set) :: ilevels    !< isentropes [K]

    ! derived variables
    !
    REAL(wp), POINTER ::      &
      &  p3d(:,:,:),          & !< 3D pressure level target field for output on p-levels
      &  z3d(:,:,:),          & !< 3D height level target field for output on z-levels
      &  i3d(:,:,:)             !< 3D theta level target field for output on isentropes

  END TYPE t_nh_pzlev_config

  !>
  !!
  TYPE(t_nh_pzlev_config), TARGET :: nh_pzlev_config(0:max_dom)


CONTAINS

  !! setup components for output on pressure/height levels and isentropes
  !!
  !! Setup of additional control variables for output on pressure/height levels 
  !! and isentropes.  
  !! These may depend on the nh_pzlev-namelist and potentially other namelists. 
  !! This routine is called, after all namelists have been read and a synoptic 
  !! consistency check has been done.
  !!
  SUBROUTINE configure_nh_pzlev( jg, nproma, npromz_c, nblks_c )
  !
    INTEGER, INTENT(IN) :: jg           !< patch 
    INTEGER, INTENT(IN) :: nproma
    INTEGER, INTENT(IN) :: npromz_c
    INTEGER, INTENT(IN) :: nblks_c      !< number of blocks

    ! Local variables
    INTEGER :: ist
    INTEGER :: nlen
    INTEGER :: z_nplev, z_nzlev, z_nilev
    INTEGER :: jb, jc, jk           ! loop indices
    !-----------------------------------------------------------------------

    z_nplev = nh_pzlev_config(jg)%plevels%nvalues
    z_nzlev = nh_pzlev_config(jg)%zlevels%nvalues
    z_nilev = nh_pzlev_config(jg)%ilevels%nvalues

    ! do status output
    !
    IF (((z_nplev > 0) .OR. (z_nzlev > 0) .OR. (z_nilev > 0)) .AND. (my_process_is_stdio())) THEN
      WRITE (0,'(a)')      " "
      WRITE (0,'(a,i0)') " Output on pressure/height levels and/or isentropes: domain ", jg
      IF (z_nplev > 0) THEN
        WRITE (0,'(a)')      " selected pressure levels: "
        WRITE (0,'(5(f10.2,","))')  nh_pzlev_config(jg)%plevels%values(1:z_nplev)
        WRITE (0,'(a)')      " "
      END IF
      IF (z_nzlev > 0) THEN
        WRITE (0,'(a)')      " selected height levels: "
        WRITE (0,'(5(f10.2,","))')  nh_pzlev_config(jg)%zlevels%values(1:z_nzlev)
        WRITE (0,'(a)')      " "
      END IF
      IF (z_nilev > 0) THEN
        WRITE (0,'(a)')      " selected isentropic levels: "
        WRITE (0,'(5(f10.2,","))')  nh_pzlev_config(jg)%ilevels%values(1:z_nilev)
        WRITE (0,'(a)')      " "
      END IF
    END IF

    ! allocate 3D pressure and z-level fields
    ALLOCATE(nh_pzlev_config(jg)%p3d(nproma,z_nplev,nblks_c),          &
      &      nh_pzlev_config(jg)%z3d(nproma,z_nzlev,nblks_c),          &
      &      nh_pzlev_config(jg)%i3d(nproma,z_nilev,nblks_c), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ( 'mo_nh_pzlev_nml: configure_nh_pzlev',       &
        &      'allocation of p3d, z3d, i3d failed' )
    ENDIF
    !$ACC ENTER DATA ASYNC(1) COPYIN(nh_pzlev_config(jg:jg))
    ! The %values components are allocated on CPU via collect_requested_ipz_levels
    !$ACC ENTER DATA ASYNC(1) &
    !$ACC   COPYIN(nh_pzlev_config(jg)%plevels%values) &
    !$ACC   COPYIN(nh_pzlev_config(jg)%zlevels%values) &
    !$ACC   COPYIN(nh_pzlev_config(jg)%ilevels%values) &
    !$ACC   CREATE(nh_pzlev_config(jg)%p3d, nh_pzlev_config(jg)%z3d, nh_pzlev_config(jg)%i3d)

    ! Fill z3d field of pressure-level data and pressure field of 
    ! height-level data
    !

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,nlen)
    ! DEFAULT(PRESENT) may check for whole nh_pzlev_config, however, as it is
    ! allocated (COPYIN) only one element at the time check explicitly only
    ! the current element
    !$ACC PARALLEL PRESENT(nh_pzlev_config(jg:jg)) DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP SEQ
    DO jb = 1,nblks_c

      IF (jb /= nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = npromz_c
      ENDIF

      !$ACC LOOP SEQ
      DO jk = 1, z_nplev
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = 1, nlen
          nh_pzlev_config(jg)%p3d(jc,jk,jb) = nh_pzlev_config(jg)%plevels%values(jk)
        ENDDO
      ENDDO

      !$ACC LOOP SEQ
      DO jk = 1, z_nzlev
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = 1, nlen
          nh_pzlev_config(jg)%z3d(jc,jk,jb) = nh_pzlev_config(jg)%zlevels%values(jk)
        ENDDO
      ENDDO

      !$ACC LOOP SEQ
      DO jk = 1, z_nilev
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = 1, nlen
          nh_pzlev_config(jg)%i3d(jc,jk,jb) = nh_pzlev_config(jg)%ilevels%values(jk)
        ENDDO
      ENDDO

    ENDDO
    !$ACC END PARALLEL
!$OMP END DO
!$OMP END PARALLEL
    
  END SUBROUTINE configure_nh_pzlev


  SUBROUTINE deallocate_nh_pzlev( jg )
    INTEGER, INTENT(IN) :: jg
    INTEGER :: ist
    CHARACTER(LEN=*), PARAMETER :: routine = "mo_nh_pzlev_nml: deallocate_nh_pzlev"
    ! deallocate 3D pressure, z-, and i-level fields

    !$ACC WAIT

    IF (ALLOCATED(nh_pzlev_config(jg)%plevels%values)) THEN
      !$ACC EXIT DATA DELETE(nh_pzlev_config(jg)%plevels%values)
      CALL deallocate_set(nh_pzlev_config(jg)%plevels)
    ENDIF

    IF (ALLOCATED(nh_pzlev_config(jg)%zlevels%values)) THEN
      !$ACC EXIT DATA DELETE(nh_pzlev_config(jg)%zlevels%values)
      CALL deallocate_set(nh_pzlev_config(jg)%zlevels)
    ENDIF

    IF (ALLOCATED(nh_pzlev_config(jg)%ilevels%values)) THEN
      !$ACC EXIT DATA DELETE(nh_pzlev_config(jg)%ilevels%values)
      CALL deallocate_set(nh_pzlev_config(jg)%ilevels)
    ENDIF

    !$ACC EXIT DATA DELETE(nh_pzlev_config(jg)%p3d, nh_pzlev_config(jg)%z3d, nh_pzlev_config(jg)%i3d)
    !$ACC EXIT DATA DELETE(nh_pzlev_config(jg:jg))
    DEALLOCATE( &
      nh_pzlev_config(jg)%p3d, &
      nh_pzlev_config(jg)%z3d, &
      nh_pzlev_config(jg)%i3d, &
      STAT=ist )
    IF (ist /= SUCCESS) CALL finish( routine, 'deallocation of p3d, z3d, i3d failed' )
  END SUBROUTINE deallocate_nh_pzlev

END MODULE mo_nh_pzlev_config
