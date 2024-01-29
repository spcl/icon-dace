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
MODULE mo_reader_cams

  USE mo_kind,                    ONLY: dp, wp, i8
  USE mo_parallel_config,         ONLY: get_nproma
  USE mo_exception,               ONLY: finish
  USE mo_reader_abstract,         ONLY: t_abstract_reader
  USE mo_impl_constants,          ONLY: n_camsaermr
  USE mo_io_units,                ONLY: FILENAME_MAX
  USE mo_model_domain,            ONLY: t_patch
  USE mo_netcdf_errhandler,       ONLY: nf
  USE mo_netcdf
  USE mtime,                      ONLY: julianday, getJulianDayFromDatetime,         &
                                    &   datetime, newdatetime, deallocatedatetime,   &
                                    &   OPERATOR(+), ASSIGNMENT(=),                  &
                                    &   no_of_ms_in_a_day, no_of_ms_in_a_hour,       &
                                    &   no_of_ms_in_a_minute, no_of_ms_in_a_second
  USE mo_mpi,                     ONLY: my_process_is_mpi_workroot, process_mpi_root_id, p_comm_work, p_bcast
  USE mo_read_netcdf_distributed, ONLY: distrib_nf_open, distrib_read, distrib_nf_close, idx_blk_time
  USE mo_fortran_tools,           ONLY: t_ptr_4d

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: t_cams_reader

  TYPE, EXTENDS(t_abstract_reader) :: t_cams_reader

    TYPE(t_patch), POINTER      :: p_patch => NULL()
    CHARACTER(len=NF_MAX_NAME)  :: varnames(n_camsaermr+1)
    CHARACTER(len=FILENAME_MAX) :: filename
    INTEGER                     :: fileid, dist_fileid,nlev_cams
    LOGICAL                     :: lopened = .FALSE.

  CONTAINS

    PROCEDURE :: init            => cams_init_reader
    PROCEDURE :: get_one_timelev => cams_get_one_timelevel
    PROCEDURE :: get_times       => cams_get_times
    PROCEDURE :: deinit          => cams_deinit_reader

    PROCEDURE :: get_nblks       => cams_get_nblks
    PROCEDURE :: get_npromz      => cams_get_npromz

  END TYPE t_cams_reader

  CHARACTER(len=*), PARAMETER :: modname = 'mo_reader_cams'

CONTAINS

  SUBROUTINE cams_init_reader(this, p_patch, filename)
    CLASS(t_cams_reader),    INTENT(inout) :: this
    TYPE(t_patch),      TARGET, INTENT(in   ) :: p_patch
    CHARACTER(len=*),           INTENT(in   ) :: filename

    CHARACTER(len=*), PARAMETER :: routine = 'cams_init_reader'

    this%filename = TRIM(filename)

    this%varnames(1)  = "Sea_Salt_bin1"              ! layer-integrated mass (kg/m2)
    this%varnames(2)  = "Sea_Salt_bin2"              ! layer-integrated mass (kg/m2)
    this%varnames(3)  = "Sea_Salt_bin3"              ! layer-integrated mass (kg/m2)
    this%varnames(4)  = "Mineral_Dust_bin1"          ! layer-integrated mass (kg/m2)
    this%varnames(5)  = "Mineral_Dust_bin2"          ! layer-integrated mass (kg/m2)
    this%varnames(6)  = "Mineral_Dust_bin3"          ! layer-integrated mass (kg/m2)
    this%varnames(7)  = "Organic_Matter_hydrophilic" ! layer-integrated mass (kg/m2)
    this%varnames(8)  = "Organic_Matter_hydrophobic" ! layer-integrated mass (kg/m2)
    this%varnames(9)  = "Black_Carbon_hydrophilic"   ! layer-integrated mass (kg/m2)
    this%varnames(10) = "Black_Carbon_hydrophobic"   ! layer-integrated mass (kg/m2)
    this%varnames(11) = "Sulfates"                   ! layer-integrated mass (kg/m2)
    this%varnames(12) = "half_level_pressure"        ! Pressure at base of layer (Pa)

    this%nlev_cams = 60
    this%p_patch => p_patch

    IF (.NOT. this%lopened) THEN
      IF (my_process_is_mpi_workroot()) THEN
        CALL nf(nf_open(this%filename, nf_nowrite, this%fileid), routine)
      ENDIF
      this%dist_fileid = distrib_nf_open(TRIM(this%filename))
      this%lopened = .TRUE.
    ENDIF

  END SUBROUTINE cams_init_reader

  SUBROUTINE cams_get_times (this, times)

    CLASS(t_cams_reader), INTENT(inout) :: &
      &  this
    TYPE(julianday), ALLOCATABLE, INTENT(out) :: &
      &  times(:)

    ! Routine is not used for climatology (implementation is still necessary)

  END SUBROUTINE cams_get_times


  SUBROUTINE cams_get_one_timelevel(this, timelevel, varname, dat)
    CLASS(t_cams_reader), INTENT(inout)  :: this
    INTEGER, INTENT(in   )               :: timelevel
    CHARACTER(len=*), INTENT(in   )      :: varname
    REAL(wp), ALLOCATABLE, INTENT(inout) :: dat(:,:,:,:)
    REAL(wp), ALLOCATABLE, TARGET        :: temp(:,:,:,:)
    TYPE(t_ptr_4d)                       :: tmp(1)
    INTEGER                              :: var_dimlen(3),var_start(3), var_end(3),jt


    ALLOCATE(temp(get_nproma(), this%nlev_cams, this%p_patch%nblks_c, n_camsaermr+1))
    IF (ALLOCATED(dat)) DEALLOCATE(dat)
    ALLOCATE( dat(get_nproma(), this%nlev_cams, this%p_patch%nblks_c, n_camsaermr+1))

      var_dimlen(2) = SIZE(temp, 2) ! number of vertical levels
      var_dimlen(3) = 1             ! number of time steps = 1
      var_start(:)  = (/1, 1, 1/)
      var_end(:)    = var_dimlen(:)
      var_start(3)  = timelevel
      var_end(3)    = timelevel

    temp(:,:,:,:) = -1.0_wp

    IF (.NOT. this%lopened) &
      CALL finish(modname, 'CAMS climatology file not open!')
    IF (TRIM(varname) /= '') &
      CALL finish(modname, 'CAMS climatology: Only bulk reading of all variables implemented!')

    tmp(1)%p => temp

    DO jt = 1, n_camsaermr+1

       CALL distrib_read(this%dist_fileid, this%varnames(jt), tmp, &
         & (/this%p_patch%cells%dist_io_data/), edim=var_dimlen(2:3), dimo=idx_blk_time, &
         & start_ext_dim=var_start(2:3), end_ext_dim=var_end(2:3))

         dat(:,:,:,jt) = temp(:,:,:,1)
    ENDDO

    DEALLOCATE(temp)
    
  END SUBROUTINE cams_get_one_timelevel

  FUNCTION cams_get_nblks (this) RESULT(nblks)
    CLASS(t_cams_reader), INTENT(in   ) :: this
    INTEGER                                :: nblks
    nblks = this%p_patch%nblks_c
  END FUNCTION cams_get_nblks

  FUNCTION cams_get_npromz (this) RESULT(npromz)
    CLASS(t_cams_reader), INTENT(in   ) :: this
    INTEGER                                :: npromz
    npromz = this%p_patch%npromz_c
  END FUNCTION cams_get_npromz

  SUBROUTINE cams_deinit_reader(this)
    CLASS(t_cams_reader), INTENT(inout) :: this
    CHARACTER(len=*), PARAMETER :: routine = 'cams_deinit_reader'
    IF (ASSOCIATED(this%p_patch)) NULLIFY(this%p_patch)
    IF (this%lopened) THEN
      IF (my_process_is_mpi_workroot()) THEN
        CALL nf(nf_close(this%fileid), routine)
      END IF
      CALL distrib_nf_close(this%dist_fileid)
    END IF
  END SUBROUTINE cams_deinit_reader

END MODULE mo_reader_cams
