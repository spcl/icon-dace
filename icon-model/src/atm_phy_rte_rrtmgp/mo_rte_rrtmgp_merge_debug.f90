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

MODULE mo_rte_rrtmgp_merge_debug

  USE mo_kind, ONLY: wp
  USE mo_netcdf

  IMPLICIT NONE

  PRIVATE

  INTEGER, PARAMETER :: ndim = 5, nvar_all = 47
  
  PUBLIC :: ndim, nvar_all

  TYPE :: Tdimension
    CHARACTER(len=32) :: name
    INTEGER :: length
  END TYPE
  TYPE(Tdimension) :: dims(ndim) = (/&
    Tdimension('kbdim', -1), &
    Tdimension('klev', -1), &
    Tdimension('klevp1', -1), &
    Tdimension('nseed', -1), &
    Tdimension('nrec', -1)/)

  INTEGER, PARAMETER :: kbdim_id = 1, klev_id = 2, klevp1_id = 3, &
    nseed_id = 4, nrec_id = 5

  TYPE :: TV !TVariable
    CHARACTER(len=32) :: name
    INTEGER :: dims(4), nf_type, ndim
    LOGICAL :: always_active
  END TYPE

  TYPE(TV), SAVE :: vars(nvar_all) = (/ &
!In write_record_interface_aes
    TV('kproma', (/5,-1,-1,-1/), nf90_int, 1, .true.), &
    TV('cosmu0', (/1,5,-1,-1/), nf90_double, 2, .true.), &
    TV('day_frc', (/1,5,-1,-1/), nf90_double, 2, .true.), &
    TV('rnseeds1', (/1,4,5,-1/), nf90_int, 3, .true.), &
    TV('rnseeds2', (/1,4,5,-1/), nf90_int, 3, .true.), &
    TV('a_ndif', (/1,5,-1,-1/), nf90_double, 2, .true.), &
    TV('a_ndir', (/1,5,-1,-1/), nf90_double, 2, .true.), &
    TV('a_vdif', (/1,5,-1,-1/), nf90_double, 2, .true.), &
    TV('a_vdir', (/1,5,-1,-1/), nf90_double, 2, .true.), &
    TV('zf', (/1,2,5,-1/), nf90_double, 3, .true.), &
    TV('zh', (/1,2,5,-1/), nf90_double, 3, .true.), &
    TV('dz', (/1,2,5,-1/), nf90_double, 3, .true.), &
    TV('pp_fl_in', (/1,2,5,-1/), nf90_double, 3, .true.), &
    TV('pp_hl_in', (/1,3,5,-1/), nf90_double, 3, .true.), &
    TV('tk_fl_in', (/1,2,5,-1/), nf90_double, 3, .true.), &
    TV('tk_hl_in', (/1,3,5,-1/), nf90_double, 3, .true.), &
    TV('pp_fl', (/1,2,5,-1/), nf90_double, 3, .true.), &
    TV('pp_hl', (/1,3,5,-1/), nf90_double, 3, .true.), &
    TV('tk_fl', (/1,2,5,-1/), nf90_double, 3, .true.), &
    TV('tk_hl', (/1,3,5,-1/), nf90_double, 3, .true.), &
    TV('tk_sfc', (/1,5,-1,-1/), nf90_double, 2, .true.), &
    TV('xm_dry', (/1,2,5,-1/), nf90_double, 3, .true.), &
    TV('xm_vap', (/1,2,5,-1/), nf90_double, 3, .true.), &
    TV('xm_co2', (/1,2,5,-1/), nf90_double, 3, .true.), &
    TV('xm_ch4', (/1,2,5,-1/), nf90_double, 3, .true.), &
    TV('xm_o2', (/1,2,5,-1/), nf90_double, 3, .true.), &
    TV('xm_o3', (/1,2,5,-1/), nf90_double, 3, .true.), &
    TV('xm_n2o', (/1,2,5,-1/), nf90_double, 3, .true.), &
    TV('cdnc', (/1,2,5,-1/), nf90_double, 3, .true.), &
    TV('cld_frc_ext', (/1,2,5,-1/), nf90_double, 3, .true.), &
    TV('f_ldcs', (/1,3,5,-1/), nf90_double, 3, .true.), &
    TV('f_lucs', (/1,3,5,-1/), nf90_double, 3, .true.), &
    TV('f_sdcs', (/1,3,5,-1/), nf90_double, 3, .true.), &
    TV('f_sucs', (/1,3,5,-1/), nf90_double, 3, .true.), &
    TV('f_ld', (/1,3,5,-1/), nf90_double, 3, .true.), &
    TV('f_lu', (/1,3,5,-1/), nf90_double, 3, .true.), &
    TV('f_sd', (/1,3,5,-1/), nf90_double, 3, .true.), &
    TV('f_su', (/1,3,5,-1/), nf90_double, 3, .true.), &
    TV('vds_dir', (/1,5,-1,-1/), nf90_double, 2, .true.), &
    TV('pds_dir', (/1,5,-1,-1/), nf90_double, 2, .true.), &
    TV('nds_dir', (/1,5,-1,-1/), nf90_double, 2, .true.), &
    TV('vds_dif', (/1,5,-1,-1/), nf90_double, 2, .true.), &
    TV('pds_dif', (/1,5,-1,-1/), nf90_double, 2, .true.), &
    TV('nds_dif', (/1,5,-1,-1/), nf90_double, 2, .true.), &
    TV('vus', (/1,5,-1,-1/), nf90_double, 2, .true.), &
    TV('pus', (/1,5,-1,-1/), nf90_double, 2, .true.), &
    TV('nus', (/1,5,-1,-1/), nf90_double, 2, .true.)/)

  TYPE :: TDumpDictionary
    INTEGER :: nvar_active
    LOGICAL :: is_active(nvar_all)
    INTEGER :: nf_idx(nvar_all)
  END TYPE

  INTERFACE crash
    MODULE PROCEDURE crash_nothing
    MODULE PROCEDURE crash_msg
    MODULE PROCEDURE crash_code
  END INTERFACE crash

  INTEGER :: nf_write_id, writing, dump_index
  LOGICAL :: dump_finished, is_omp_radiation
#ifdef _OPENMP
!$OMP THREADPRIVATE(dump_index)
#endif
  INTEGER, SAVE :: dump_offset = 0
  INTEGER, PARAMETER :: ngas_old = 6

  TYPE(TDumpDictionary) :: dump_list

  INTEGER, PARAMETER :: dump_max = 16384
  INTEGER :: current_column

  PUBLIC :: dump_offset, dump_index, dump_finished, dump_max, is_omp_radiation
  PUBLIC :: open_write_icon
  PUBLIC :: close_write, write_record_interface_aes

  PUBLIC :: init_dump_dictionary, dim_id, var_id, vars, dims, TDumpDictionary

  PUBLIC :: init_vars_ndim
  PUBLIC :: current_column

  PUBLIC :: crash


CONTAINS

  SUBROUTINE init_vars_ndim
    INTEGER :: i
    DO i = 1, nvar_all
      vars(i)%ndim = COUNT(vars(i)%dims > 0)
    END DO
  END SUBROUTINE init_vars_ndim

  SUBROUTINE init_dump_dictionary(dic)
    TYPE(TDumpDictionary), INTENT(INOUT) :: dic
    INTEGER :: id
    dic%nvar_active = 0
    DO id = 1, nvar_all
      dic%is_active(id) = vars(id)%always_active
    ENDDO
    dic%nf_idx = -1
  END SUBROUTINE init_dump_dictionary

  INTEGER FUNCTION var_id(name)
    CHARACTER(len=*) :: name
    INTEGER :: i

    var_id = 0
    DO i = 1, nvar_all
      IF (TRIM(name) == TRIM(vars(i)%name)) THEN
        var_id = i
        EXIT
      ENDIF
    ENDDO
  END FUNCTION

  INTEGER FUNCTION dim_id(name)
    CHARACTER(len=*) :: name
    INTEGER :: i

    dim_id = 0
    DO i = 1, nvar_all
      IF (TRIM(name) == TRIM(dims(i)%name)) THEN
        dim_id = i
        EXIT
      ENDIF
    ENDDO
  END FUNCTION

  SUBROUTINE open_write_icon(kbdim, klev)
    USE mo_mpi, ONLY: get_my_mpi_all_id
    USE mo_time_config, ONLY: restart_ini_datetime_string
    USE mo_master_config, ONLY: isRestart
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: kbdim, klev
    
    INTEGER :: p_pe
    CHARACTER(len=256) :: filename

    IF (writing == 1 .or. dump_finished) THEN
      RETURN 
    END IF
    writing = 1

    CALL init_dump_dictionary(dump_list)
    p_pe = get_my_mpi_all_id()
    IF (isRestart()) THEN
      WRITE (filename,'(a,i3.3,a,a,a)') 'psdump_', p_pe, &
        '_', TRIM(restart_ini_datetime_string), '.nc'
    ELSE
      WRITE (filename,'(a,i3.3,a)') 'psdump_', p_pe, '.nc'
    ENDIF
    CALL open_write_internal(filename, kbdim, klev)
  END SUBROUTINE open_write_icon

  SUBROUTINE open_write_internal(name, kbdim, klev)
    IMPLICIT NONE

    CHARACTER(len=*), INTENT(IN) :: name
    INTEGER, INTENT(IN) :: kbdim, klev
    
    INTEGER :: i, ret, dim_val(ndim), idummy

    ret = nf90_create(TRIM(name), nf90_clobber, nf_write_id)
    IF (ret /= nf90_noerr) CALL crash(ret, __LINE__)

    dim_val = (/kbdim, klev, klev+1, 4, nf90_unlimited/)
    DO i = 1,ndim
      dims(i)%length = dim_val(i)
      ret = nf90_def_dim(nf_write_id, TRIM(dims(i)%name), dim_val(i), idummy)
      IF (ret /= nf90_noerr) CALL crash(ret, __LINE__)
    END DO

    DO i = 1, nvar_all
      IF (dump_list%is_active(i)) THEN
        ret = nf90_def_var(nf_write_id, TRIM(vars(i)%name), &
          vars(i)%nf_type, vars(i)%dims, dump_list%nf_idx(i))
        IF (ret /= nf90_noerr) CALL crash(ret, __LINE__)
      ENDIF
    END DO

    ret = nf90_enddef(nf_write_id)
    if (ret /= nf90_noerr) CALL crash(ret, __LINE__)

  END SUBROUTINE open_write_internal

  SUBROUTINE write_double(var_name, data, pos)
    CHARACTER(len=*) :: var_name
    INTEGER, INTENT(IN) :: pos
    REAL(wp), INTENT(IN) :: data(*)
    INTEGER :: i, ret, var_idx, start(ndim), count(ndim)

    var_idx = var_id(var_name)
    IF (var_idx == 0) THEN
      CALL crash('Variable ' // TRIM(var_name) // ' not in dump dictionary', &
        __LINE__)
    ENDIF

    IF (.not. dump_list%is_active(var_idx)) RETURN
    IF (vars(var_idx)%nf_type /= nf90_double) THEN
      CALL crash('Wrong type in call', __LINE__)
    END IF
    start = 1
    DO i = 1,vars(var_idx)%ndim-1
      count(i) = dims(vars(var_idx)%dims(i))%length
    END DO
    count(i) = 1
    start(i) = pos
    ret = nf90_put_var(nf_write_id, &
      dump_list%nf_idx(var_idx), data, start, count)
    if (ret /= nf90_noerr) CALL crash(ret, __LINE__)
  END SUBROUTINE write_double

  SUBROUTINE write_int(var_name, data, pos)
    CHARACTER(len=*) :: var_name
    INTEGER, INTENT(IN) :: pos
    INTEGER, INTENT(IN) :: data(*)
    INTEGER :: i, ret, var_idx, start(ndim), count(ndim)

    var_idx = var_id(var_name)
    IF (var_idx == 0) THEN
      CALL crash('Variable ' // TRIM(var_name) // ' not in dump dictionary', &
        __LINE__)
    ENDIF

    IF (.not. dump_list%is_active(var_idx)) RETURN
    IF (vars(var_idx)%nf_type /= nf90_int) THEN
      CALL crash('Wrong type in call', __LINE__)
    END IF
    start = 1
    DO i = 1,vars(var_idx)%ndim-1
      count(i) = dims(vars(var_idx)%dims(i))%length
    END DO
    count(i) = 1
    start(i) = pos
    ret = nf90_put_var(nf_write_id, &
      dump_list%nf_idx(var_idx), data, start, count)
    if (ret /= nf90_noerr) CALL crash(ret, __LINE__)
  END SUBROUTINE write_int

  SUBROUTINE write_record_interface_aes(kproma, cosmu0, day_frc, &
    rnseeds1, rnseeds2, &
    a_vdir, a_ndir, a_vdif, a_ndif, tk_sfc, zf, zh, dz, &
    pp_fl_in, pp_hl_in, tk_fl_in, tk_hl_in, &
    pp_fl, pp_hl, tk_fl, tk_hl, &
    xm_dry, xm_vap, xm_co2, xm_ch4, xm_o2, xm_o3, xm_n2o, &
    cdnc, cld_frc_ext, &
    f_ldcs, f_lucs, f_sdcs, f_sucs, f_ld, f_lu, f_sd, f_su, &
    vds_dir, pds_dir, nds_dir, vds_dif, pds_dif, nds_dif, vus, pus, nus)

    INTEGER, INTENT(IN) :: kproma, rnseeds1(:,:),  rnseeds2(:,:)

    REAL(wp), DIMENSION(:), INTENT(IN) :: cosmu0, day_frc, &
      a_vdir, a_ndir, a_vdif, a_ndif, tk_sfc, &
      vds_dir, pds_dir, nds_dir, vds_dif, pds_dif, nds_dif, vus, pus, nus
    
    REAL(wp), DIMENSION(:,:), INTENT(IN) :: &
      zf, zh, dz, pp_fl_in, pp_hl_in, tk_fl_in, tk_hl_in, &
      pp_fl, pp_hl, tk_fl, tk_hl, &
      xm_dry, xm_vap, xm_co2, xm_ch4, xm_o2, xm_o3, xm_n2o, &
      cdnc, cld_frc_ext, &
      f_ldcs, f_lucs, f_sdcs, f_sucs, f_ld, f_lu, f_sd, f_su

    IF (dump_finished .or. dump_index > dump_max) THEN
      RETURN
    END IF

    CALL write_int('kproma', [kproma], dump_index)
    CALL write_double('cosmu0', cosmu0, dump_index)
    CALL write_double('day_frc', day_frc, dump_index)
    CALL write_int('rnseeds1', rnseeds1, dump_index)
    CALL write_int('rnseeds2', rnseeds2, dump_index)
    CALL write_double('a_ndif', a_ndif, dump_index)
    CALL write_double('a_ndir', a_ndir, dump_index)
    CALL write_double('a_vdif', a_vdif, dump_index)
    CALL write_double('a_vdir', a_vdir, dump_index)
    CALL write_double('zf', zf, dump_index)
    CALL write_double('zh', zh, dump_index)
    CALL write_double('dz', dz, dump_index)
    CALL write_double('pp_fl_in', pp_fl_in, dump_index)
    CALL write_double('pp_hl_in', pp_hl_in, dump_index)
    CALL write_double('tk_fl_in', tk_fl_in, dump_index)
    CALL write_double('tk_hl_in', tk_hl_in, dump_index)
    CALL write_double('pp_fl', pp_fl, dump_index)
    CALL write_double('pp_hl', pp_hl, dump_index)
    CALL write_double('tk_fl', tk_fl, dump_index)
    CALL write_double('tk_hl', tk_hl, dump_index)
    CALL write_double('tk_sfc', tk_sfc, dump_index)
    CALL write_double('xm_dry', xm_dry, dump_index)
    CALL write_double('xm_vap', xm_vap, dump_index)
    CALL write_double('xm_co2', xm_co2, dump_index)
    CALL write_double('xm_ch4', xm_ch4, dump_index)
    CALL write_double('xm_o2', xm_o2, dump_index)
    CALL write_double('xm_o3', xm_o3, dump_index)
    CALL write_double('xm_n2o', xm_n2o, dump_index)
    CALL write_double('cdnc', cdnc, dump_index)
    CALL write_double('cld_frc_ext', cld_frc_ext, dump_index)
    CALL write_double('f_ldcs', f_ldcs, dump_index)
    CALL write_double('f_lucs', f_lucs, dump_index)
    CALL write_double('f_sdcs', f_sdcs, dump_index)
    CALL write_double('f_sucs', f_sucs, dump_index)
    CALL write_double('f_ld', f_ld, dump_index)
    CALL write_double('f_lu', f_lu, dump_index)
    CALL write_double('f_sd', f_sd, dump_index)
    CALL write_double('f_su', f_su, dump_index)
    CALL write_double('vds_dir', vds_dir, dump_index)
    CALL write_double('pds_dir', pds_dir, dump_index)
    CALL write_double('nds_dir', nds_dir, dump_index)
    CALL write_double('vds_dif', vds_dif, dump_index)
    CALL write_double('pds_dif', pds_dif, dump_index)
    CALL write_double('nds_dif', nds_dif, dump_index)
    CALL write_double('vus', vus, dump_index)
    CALL write_double('pus', pus, dump_index)
    CALL write_double('nus', nus, dump_index)

  END SUBROUTINE write_record_interface_aes

  SUBROUTINE close_write
    INTEGER :: ret
    IF (writing == 1) THEN
      ret = nf90_close(nf_write_id)
      IF (ret /= nf90_noerr) CALL crash(ret, __LINE__)
      writing = 0
    ENDIF
  END SUBROUTINE close_write

  SUBROUTINE crash_nothing(msg)
    CHARACTER(len=*), INTENT(IN) :: msg

    WRITE(*,*) msg
    STOP
  END SUBROUTINE crash_nothing

  SUBROUTINE crash_msg(msg, lineno)
    CHARACTER(len=*), INTENT(IN) :: msg
    INTEGER, INTENT(IN) :: lineno

    WRITE(*,*) msg, ' at yarff.f90:', lineno
    STOP
  END SUBROUTINE crash_msg

  SUBROUTINE crash_code(code, lineno)
    INTEGER, INTENT(IN) :: code, lineno

    WRITE(*,*) 'Exit due to netcdf error ', code, &
      ' at yarrf.f90:', lineno
    STOP
  END SUBROUTINE crash_code

END MODULE mo_rte_rrtmgp_merge_debug
