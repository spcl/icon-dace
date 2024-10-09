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

! Write and read RBF coefficients

MODULE mo_ser_rbf_coefficients

#ifdef SERIALIZE
  USE m_serialize
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_grid_config,         ONLY: n_dom, n_dom_start
  USE mo_exception,           ONLY: message,message_text,warning
  USE mo_mpi,                 ONLY: get_my_mpi_work_id

  PUBLIC :: ser_rbf_coefficients

  CONTAINS

  SUBROUTINE ser_rbf_coefficients(ptr_int_state)

  TYPE(t_int_state), INTENT(INOUT) :: ptr_int_state(n_dom_start:n_dom)

  TYPE(t_serializer)               :: ser_ref
  TYPE(t_savepoint)                :: savepoint
  CHARACTER(len=100) :: suffix,savepoint_name, serializer_name
  CHARACTER(len=1)   :: mode

  INTEGER :: jg

#ifdef SERIALIZE_CREATE_REFERENCE
  mode = 'w'
#elif SERIALIZE_READ_REFERENCE
  mode = 'r'
#endif

  WRITE(serializer_name,'(A,I2.2)') 'rbf_coefficients_rank',get_my_mpi_work_id()

  CALL fs_create_serializer( 'ser_rbf', serializer_name, mode, ser_ref )

  CALL warning('SER:'//TRIM(serializer_name),'Serialization is active!')

  DO jg=n_dom_start,n_dom

    WRITE(message_text,'(A,I2.2)') 'serialize rbf coefficients for dom: ',jg
    CALL message('mo_ser_rbf_coefficients',message_text)

    savepoint_name = 'rbf'
    WRITE(suffix,'(A,I2.2)') '_dom',jg
    CALL fs_create_savepoint(TRIM(savepoint_name)//TRIM(suffix), savepoint)
    CALL fs_add_savepoint_metainfo(savepoint, 'domain', jg)
    CALL fs_add_savepoint_metainfo(savepoint, 'rank', get_my_mpi_work_id())

#ifdef SERIALIZE_CREATE_REFERENCE
    CALL fs_write_field( ser_ref, savepoint, 'rbf_vec_idx'//TRIM(suffix), ptr_int_state(jg)%rbf_vec_idx_c(:,:,:) )
    CALL fs_write_field( ser_ref, savepoint, 'rbf_vec_blk'//TRIM(suffix), ptr_int_state(jg)%rbf_vec_blk_c(:,:,:) )
    CALL fs_write_field( ser_ref, savepoint, 'rbf_vec_stencil'//TRIM(suffix), ptr_int_state(jg)%rbf_vec_stencil_c(:,:) )
    CALL fs_write_field( ser_ref, savepoint, 'rbf_vec_coeff'//TRIM(suffix), ptr_int_state(jg)%rbf_vec_coeff_c(:,:,:,:) )
    CALL fs_write_field( ser_ref, savepoint, 'rbf_c2grad_idx'//TRIM(suffix), ptr_int_state(jg)%rbf_c2grad_idx(:,:,:) )
    CALL fs_write_field( ser_ref, savepoint, 'rbf_c2grad_blk'//TRIM(suffix), ptr_int_state(jg)%rbf_c2grad_blk(:,:,:) )
    CALL fs_write_field( ser_ref, savepoint, 'rbf_c2grad_coeff'//TRIM(suffix), ptr_int_state(jg)%rbf_c2grad_coeff(:,:,:,:) )
    CALL fs_write_field( ser_ref, savepoint, 'rbf_vec_idx_v'//TRIM(suffix), ptr_int_state(jg)%rbf_vec_idx_v(:,:,:) )
    CALL fs_write_field( ser_ref, savepoint, 'rbf_vec_blk_v'//TRIM(suffix), ptr_int_state(jg)%rbf_vec_blk_v(:,:,:) )
    CALL fs_write_field( ser_ref, savepoint, 'rbf_vec_stencil_v'//TRIM(suffix), ptr_int_state(jg)%rbf_vec_stencil_v(:,:) )
    CALL fs_write_field( ser_ref, savepoint, 'rbf_vec_coeff_v'//TRIM(suffix), ptr_int_state(jg)%rbf_vec_coeff_v(:,:,:,:) )
    CALL fs_write_field( ser_ref, savepoint, 'rbf_vec_idx_e'//TRIM(suffix), ptr_int_state(jg)%rbf_vec_idx_e(:,:,:) )
    CALL fs_write_field( ser_ref, savepoint, 'rbf_vec_blk_e'//TRIM(suffix), ptr_int_state(jg)%rbf_vec_blk_e(:,:,:) )
    CALL fs_write_field( ser_ref, savepoint, 'rbf_vec_stencil_e'//TRIM(suffix), ptr_int_state(jg)%rbf_vec_stencil_e(:,:) )
    CALL fs_write_field( ser_ref, savepoint, 'rbf_vec_coeff_e'//TRIM(suffix), ptr_int_state(jg)%rbf_vec_coeff_e(:,:,:) )
#elif SERIALIZE_READ_REFERENCE
    CALL fs_read_field( ser_ref, savepoint, 'rbf_vec_idx'//TRIM(suffix), ptr_int_state(jg)%rbf_vec_idx_c(:,:,:) )
    CALL fs_read_field( ser_ref, savepoint, 'rbf_vec_blk'//TRIM(suffix), ptr_int_state(jg)%rbf_vec_blk_c(:,:,:) )
    CALL fs_read_field( ser_ref, savepoint, 'rbf_vec_stencil'//TRIM(suffix), ptr_int_state(jg)%rbf_vec_stencil_c(:,:) )
    CALL fs_read_field( ser_ref, savepoint, 'rbf_vec_coeff'//TRIM(suffix), ptr_int_state(jg)%rbf_vec_coeff_c(:,:,:,:) )
    CALL fs_read_field( ser_ref, savepoint, 'rbf_c2grad_idx'//TRIM(suffix), ptr_int_state(jg)%rbf_c2grad_idx(:,:,:) )
    CALL fs_read_field( ser_ref, savepoint, 'rbf_c2grad_blk'//TRIM(suffix), ptr_int_state(jg)%rbf_c2grad_blk(:,:,:) )
    CALL fs_read_field( ser_ref, savepoint, 'rbf_c2grad_coeff'//TRIM(suffix), ptr_int_state(jg)%rbf_c2grad_coeff(:,:,:,:) )
    CALL fs_read_field( ser_ref, savepoint, 'rbf_vec_idx_v'//TRIM(suffix), ptr_int_state(jg)%rbf_vec_idx_v(:,:,:) )
    CALL fs_read_field( ser_ref, savepoint, 'rbf_vec_blk_v'//TRIM(suffix), ptr_int_state(jg)%rbf_vec_blk_v(:,:,:) )
    CALL fs_read_field( ser_ref, savepoint, 'rbf_vec_stencil_v'//TRIM(suffix), ptr_int_state(jg)%rbf_vec_stencil_v(:,:) )
    CALL fs_read_field( ser_ref, savepoint, 'rbf_vec_coeff_v'//TRIM(suffix), ptr_int_state(jg)%rbf_vec_coeff_v(:,:,:,:) )
    CALL fs_read_field( ser_ref, savepoint, 'rbf_vec_idx_e'//TRIM(suffix), ptr_int_state(jg)%rbf_vec_idx_e(:,:,:) )
    CALL fs_read_field( ser_ref, savepoint, 'rbf_vec_blk_e'//TRIM(suffix), ptr_int_state(jg)%rbf_vec_blk_e(:,:,:) )
    CALL fs_read_field( ser_ref, savepoint, 'rbf_vec_stencil_e'//TRIM(suffix), ptr_int_state(jg)%rbf_vec_stencil_e(:,:) )
    CALL fs_read_field( ser_ref, savepoint, 'rbf_vec_coeff_e'//TRIM(suffix), ptr_int_state(jg)%rbf_vec_coeff_e(:,:,:) )
#endif

    CALL fs_destroy_savepoint( savepoint )
  ENDDO

  CALL fs_destroy_serializer( ser_ref )

  END SUBROUTINE ser_rbf_coefficients

#endif

END MODULE mo_ser_rbf_coefficients
