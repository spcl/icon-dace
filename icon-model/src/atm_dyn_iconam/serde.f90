MODULE serde
  IMPLICIT NONE
  INTERFACE serialize
    MODULE PROCEDURE :: W_string
    MODULE PROCEDURE W_t_grid_domain_decomp_info, W_t_int_state, W_t_grid_cells, W_t_grid_edges, W_t_grid_vertices, W_t_patch, W_t_nh_prog, W_t_nh_diag, W_t_nh_metrics, W_logical_R_1, W_integer__1_R_1, W_integer__2_R_1, W_integer__4_R_1, W_integer__8_R_1, W_real__4_R_1, W_real__8_R_1, W_logical_R_2, W_integer__1_R_2, W_integer__2_R_2, W_integer__4_R_2, W_integer__8_R_2, W_real__4_R_2, W_real__8_R_2, W_logical_R_3, W_integer__1_R_3, W_integer__2_R_3, W_integer__4_R_3, W_integer__8_R_3, W_real__4_R_3, W_real__8_R_3, W_logical_R_4, W_integer__1_R_4, W_integer__2_R_4, W_integer__4_R_4, W_integer__8_R_4, W_real__4_R_4, W_real__8_R_4, W_logical, W_integer1, W_integer2, W_integer4, W_integer8, W_real4, W_real8
  END INTERFACE serialize
  INTEGER :: generation = 0
  CONTAINS
  SUBROUTINE tic
    generation = generation + 1
  END SUBROUTINE tic
  FUNCTION cat(prefix, asis) RESULT(path)
    CHARACTER(LEN = *), INTENT(IN) :: prefix
    CHARACTER(LEN = :), ALLOCATABLE :: path
    CHARACTER(LEN = 50) :: gen
    LOGICAL, INTENT(IN) :: asis
    IF (asis) THEN
      path = prefix
    ELSE
      WRITE(gen, '(g0)') generation
      path = prefix // '.' // TRIM(gen) // ".data"
    END IF
  END FUNCTION cat
  FUNCTION at(prefix, asis) RESULT(io)
    CHARACTER(LEN = *), INTENT(IN) :: prefix
    INTEGER :: io
    LOGICAL, OPTIONAL, INTENT(IN) :: asis
    LOGICAL :: asis_local
    asis_local = .FALSE.
    IF (PRESENT(asis)) asis_local = asis
    OPEN(NEWUNIT = io, FILE = cat(prefix, asis_local), STATUS = "replace", ACTION = "write")
  END FUNCTION at
  SUBROUTINE W_string(io, x, cleanup, nline)
    INTEGER :: io
    CHARACTER(LEN = *), INTENT(IN) :: x
    INTEGER :: i, xend
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline
    LOGICAL :: cleanup_local, nline_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    xend = LEN(x)
    DO i = 1, LEN(x)
      IF (x(i : i) == CHAR(0)) THEN
        xend = i - 1
        EXIT
      END IF
    END DO
    WRITE(io, '(A)', ADVANCE = 'no') TRIM(x(1 : xend))
    IF (nline_local) WRITE(io, '(g0)', ADVANCE = 'no') NEW_LINE('A')
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_string
  SUBROUTINE W_t_grid_domain_decomp_info(io, x, cleanup, nline)
    USE mo_decomposition_tools, ONLY: t_grid_domain_decomp_info
    INTEGER :: io
    TYPE(t_grid_domain_decomp_info), TARGET, INTENT(IN) :: x
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline
    INTEGER :: kmeta, kmeta_0, kmeta_1, kmeta_2, kmeta_3, kmeta_4, kmeta_5, kmeta_6, kmeta_7, kmeta_8, kmeta_9
    LOGICAL :: cleanup_local, nline_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    CALL serialize(io, '# owner_mask', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, ALLOCATED(x % owner_mask), cleanup = .FALSE.)
    IF (ALLOCATED(x % owner_mask)) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 2, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 2
        CALL serialize(io, SIZE(x % owner_mask, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 2
        CALL serialize(io, LBOUND(x % owner_mask, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, x % owner_mask, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    END IF
    IF (nline_local) WRITE(io, '(g0)', ADVANCE = 'no') NEW_LINE('A')
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_t_grid_domain_decomp_info
  SUBROUTINE W_t_int_state(io, x, cleanup, nline)
    USE mo_intp_data_strc, ONLY: t_int_state
    INTEGER :: io
    TYPE(t_int_state), TARGET, INTENT(IN) :: x
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline
    INTEGER :: kmeta, kmeta_0, kmeta_1, kmeta_2, kmeta_3, kmeta_4, kmeta_5, kmeta_6, kmeta_7, kmeta_8, kmeta_9
    LOGICAL :: cleanup_local, nline_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    CALL serialize(io, '# c_lin_e', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, ALLOCATED(x % c_lin_e), cleanup = .FALSE.)
    IF (ALLOCATED(x % c_lin_e)) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 3, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 3
        CALL serialize(io, SIZE(x % c_lin_e, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 3
        CALL serialize(io, LBOUND(x % c_lin_e, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, x % c_lin_e, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    END IF
    CALL serialize(io, '# e_bln_c_s', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, ALLOCATED(x % e_bln_c_s), cleanup = .FALSE.)
    IF (ALLOCATED(x % e_bln_c_s)) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 3, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 3
        CALL serialize(io, SIZE(x % e_bln_c_s, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 3
        CALL serialize(io, LBOUND(x % e_bln_c_s, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, x % e_bln_c_s, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    END IF
    CALL serialize(io, '# cells_aw_verts', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, ALLOCATED(x % cells_aw_verts), cleanup = .FALSE.)
    IF (ALLOCATED(x % cells_aw_verts)) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 3, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 3
        CALL serialize(io, SIZE(x % cells_aw_verts, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 3
        CALL serialize(io, LBOUND(x % cells_aw_verts, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, x % cells_aw_verts, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    END IF
    CALL serialize(io, '# rbf_vec_coeff_e', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, ALLOCATED(x % rbf_vec_coeff_e), cleanup = .FALSE.)
    IF (ALLOCATED(x % rbf_vec_coeff_e)) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 3, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 3
        CALL serialize(io, SIZE(x % rbf_vec_coeff_e, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 3
        CALL serialize(io, LBOUND(x % rbf_vec_coeff_e, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, x % rbf_vec_coeff_e, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    END IF
    CALL serialize(io, '# geofac_grdiv', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, ALLOCATED(x % geofac_grdiv), cleanup = .FALSE.)
    IF (ALLOCATED(x % geofac_grdiv)) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 3, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 3
        CALL serialize(io, SIZE(x % geofac_grdiv, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 3
        CALL serialize(io, LBOUND(x % geofac_grdiv, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, x % geofac_grdiv, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    END IF
    CALL serialize(io, '# geofac_rot', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, ALLOCATED(x % geofac_rot), cleanup = .FALSE.)
    IF (ALLOCATED(x % geofac_rot)) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 3, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 3
        CALL serialize(io, SIZE(x % geofac_rot, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 3
        CALL serialize(io, LBOUND(x % geofac_rot, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, x % geofac_rot, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    END IF
    CALL serialize(io, '# geofac_n2s', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, ALLOCATED(x % geofac_n2s), cleanup = .FALSE.)
    IF (ALLOCATED(x % geofac_n2s)) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 3, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 3
        CALL serialize(io, SIZE(x % geofac_n2s, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 3
        CALL serialize(io, LBOUND(x % geofac_n2s, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, x % geofac_n2s, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    END IF
    IF (nline_local) WRITE(io, '(g0)', ADVANCE = 'no') NEW_LINE('A')
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_t_int_state
  SUBROUTINE W_t_grid_cells(io, x, cleanup, nline)
    USE mo_model_domain, ONLY: t_grid_cells
    INTEGER :: io
    TYPE(t_grid_cells), TARGET, INTENT(IN) :: x
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline
    INTEGER :: kmeta, kmeta_0, kmeta_1, kmeta_2, kmeta_3, kmeta_4, kmeta_5, kmeta_6, kmeta_7, kmeta_8, kmeta_9
    LOGICAL :: cleanup_local, nline_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    CALL serialize(io, '# neighbor_idx', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, ALLOCATED(x % neighbor_idx), cleanup = .FALSE.)
    IF (ALLOCATED(x % neighbor_idx)) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 3, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 3
        CALL serialize(io, SIZE(x % neighbor_idx, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 3
        CALL serialize(io, LBOUND(x % neighbor_idx, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, x % neighbor_idx, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    END IF
    CALL serialize(io, '# neighbor_blk', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, ALLOCATED(x % neighbor_blk), cleanup = .FALSE.)
    IF (ALLOCATED(x % neighbor_blk)) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 3, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 3
        CALL serialize(io, SIZE(x % neighbor_blk, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 3
        CALL serialize(io, LBOUND(x % neighbor_blk, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, x % neighbor_blk, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    END IF
    CALL serialize(io, '# edge_idx', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, ALLOCATED(x % edge_idx), cleanup = .FALSE.)
    IF (ALLOCATED(x % edge_idx)) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 3, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 3
        CALL serialize(io, SIZE(x % edge_idx, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 3
        CALL serialize(io, LBOUND(x % edge_idx, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, x % edge_idx, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    END IF
    CALL serialize(io, '# edge_blk', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, ALLOCATED(x % edge_blk), cleanup = .FALSE.)
    IF (ALLOCATED(x % edge_blk)) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 3, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 3
        CALL serialize(io, SIZE(x % edge_blk, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 3
        CALL serialize(io, LBOUND(x % edge_blk, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, x % edge_blk, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    END IF
    CALL serialize(io, '# area', cleanup = .FALSE.)
    CALL serialize(io, '# assoc', cleanup = .FALSE.)
    CALL serialize(io, ASSOCIATED(x % area), cleanup = .FALSE.)
    IF (ASSOCIATED(x % area)) THEN
      kmeta = 0
      CALL serialize(io, "# missing", cleanup = .FALSE.)
      CALL serialize(io, (kmeta == 0), cleanup = .FALSE.)
      CALL serialize(io, x % area, cleanup = .FALSE.)
    END IF
    CALL serialize(io, '# start_index', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, ALLOCATED(x % start_index), cleanup = .FALSE.)
    IF (ALLOCATED(x % start_index)) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 1, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 1
        CALL serialize(io, SIZE(x % start_index, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 1
        CALL serialize(io, LBOUND(x % start_index, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, x % start_index, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    END IF
    CALL serialize(io, '# end_index', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, ALLOCATED(x % end_index), cleanup = .FALSE.)
    IF (ALLOCATED(x % end_index)) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 1, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 1
        CALL serialize(io, SIZE(x % end_index, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 1
        CALL serialize(io, LBOUND(x % end_index, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, x % end_index, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    END IF
    CALL serialize(io, '# start_block', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, ALLOCATED(x % start_block), cleanup = .FALSE.)
    IF (ALLOCATED(x % start_block)) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 1, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 1
        CALL serialize(io, SIZE(x % start_block, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 1
        CALL serialize(io, LBOUND(x % start_block, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, x % start_block, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    END IF
    CALL serialize(io, '# end_block', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, ALLOCATED(x % end_block), cleanup = .FALSE.)
    IF (ALLOCATED(x % end_block)) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 1, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 1
        CALL serialize(io, SIZE(x % end_block, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 1
        CALL serialize(io, LBOUND(x % end_block, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, x % end_block, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    END IF
    CALL serialize(io, '# decomp_info', cleanup = .FALSE.)
    CALL serialize(io, x % decomp_info, cleanup = .FALSE.)
    IF (nline_local) WRITE(io, '(g0)', ADVANCE = 'no') NEW_LINE('A')
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_t_grid_cells
  SUBROUTINE W_t_grid_edges(io, x, cleanup, nline)
    USE mo_model_domain, ONLY: t_grid_edges
    INTEGER :: io
    TYPE(t_grid_edges), TARGET, INTENT(IN) :: x
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline
    INTEGER :: kmeta, kmeta_0, kmeta_1, kmeta_2, kmeta_3, kmeta_4, kmeta_5, kmeta_6, kmeta_7, kmeta_8, kmeta_9
    LOGICAL :: cleanup_local, nline_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    CALL serialize(io, '# cell_idx', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, ALLOCATED(x % cell_idx), cleanup = .FALSE.)
    IF (ALLOCATED(x % cell_idx)) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 3, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 3
        CALL serialize(io, SIZE(x % cell_idx, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 3
        CALL serialize(io, LBOUND(x % cell_idx, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, x % cell_idx, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    END IF
    CALL serialize(io, '# cell_blk', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, ALLOCATED(x % cell_blk), cleanup = .FALSE.)
    IF (ALLOCATED(x % cell_blk)) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 3, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 3
        CALL serialize(io, SIZE(x % cell_blk, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 3
        CALL serialize(io, LBOUND(x % cell_blk, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, x % cell_blk, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    END IF
    CALL serialize(io, '# vertex_idx', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, ALLOCATED(x % vertex_idx), cleanup = .FALSE.)
    IF (ALLOCATED(x % vertex_idx)) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 3, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 3
        CALL serialize(io, SIZE(x % vertex_idx, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 3
        CALL serialize(io, LBOUND(x % vertex_idx, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, x % vertex_idx, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    END IF
    CALL serialize(io, '# vertex_blk', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, ALLOCATED(x % vertex_blk), cleanup = .FALSE.)
    IF (ALLOCATED(x % vertex_blk)) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 3, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 3
        CALL serialize(io, SIZE(x % vertex_blk, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 3
        CALL serialize(io, LBOUND(x % vertex_blk, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, x % vertex_blk, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    END IF
    CALL serialize(io, '# tangent_orientation', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, ALLOCATED(x % tangent_orientation), cleanup = .FALSE.)
    IF (ALLOCATED(x % tangent_orientation)) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 2, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 2
        CALL serialize(io, SIZE(x % tangent_orientation, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 2
        CALL serialize(io, LBOUND(x % tangent_orientation, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, x % tangent_orientation, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    END IF
    CALL serialize(io, '# quad_idx', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, ALLOCATED(x % quad_idx), cleanup = .FALSE.)
    IF (ALLOCATED(x % quad_idx)) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 3, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 3
        CALL serialize(io, SIZE(x % quad_idx, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 3
        CALL serialize(io, LBOUND(x % quad_idx, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, x % quad_idx, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    END IF
    CALL serialize(io, '# quad_blk', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, ALLOCATED(x % quad_blk), cleanup = .FALSE.)
    IF (ALLOCATED(x % quad_blk)) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 3, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 3
        CALL serialize(io, SIZE(x % quad_blk, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 3
        CALL serialize(io, LBOUND(x % quad_blk, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, x % quad_blk, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    END IF
    CALL serialize(io, '# inv_primal_edge_length', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, ALLOCATED(x % inv_primal_edge_length), cleanup = .FALSE.)
    IF (ALLOCATED(x % inv_primal_edge_length)) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 2, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 2
        CALL serialize(io, SIZE(x % inv_primal_edge_length, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 2
        CALL serialize(io, LBOUND(x % inv_primal_edge_length, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, x % inv_primal_edge_length, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    END IF
    CALL serialize(io, '# inv_dual_edge_length', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, ALLOCATED(x % inv_dual_edge_length), cleanup = .FALSE.)
    IF (ALLOCATED(x % inv_dual_edge_length)) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 2, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 2
        CALL serialize(io, SIZE(x % inv_dual_edge_length, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 2
        CALL serialize(io, LBOUND(x % inv_dual_edge_length, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, x % inv_dual_edge_length, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    END IF
    CALL serialize(io, '# area_edge', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, ALLOCATED(x % area_edge), cleanup = .FALSE.)
    IF (ALLOCATED(x % area_edge)) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 2, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 2
        CALL serialize(io, SIZE(x % area_edge, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 2
        CALL serialize(io, LBOUND(x % area_edge, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, x % area_edge, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    END IF
    CALL serialize(io, '# f_e', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, ALLOCATED(x % f_e), cleanup = .FALSE.)
    IF (ALLOCATED(x % f_e)) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 2, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 2
        CALL serialize(io, SIZE(x % f_e, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 2
        CALL serialize(io, LBOUND(x % f_e, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, x % f_e, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    END IF
    CALL serialize(io, '# fn_e', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, ALLOCATED(x % fn_e), cleanup = .FALSE.)
    IF (ALLOCATED(x % fn_e)) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 2, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 2
        CALL serialize(io, SIZE(x % fn_e, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 2
        CALL serialize(io, LBOUND(x % fn_e, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, x % fn_e, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    END IF
    CALL serialize(io, '# ft_e', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, ALLOCATED(x % ft_e), cleanup = .FALSE.)
    IF (ALLOCATED(x % ft_e)) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 2, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 2
        CALL serialize(io, SIZE(x % ft_e, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 2
        CALL serialize(io, LBOUND(x % ft_e, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, x % ft_e, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    END IF
    CALL serialize(io, '# start_index', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, ALLOCATED(x % start_index), cleanup = .FALSE.)
    IF (ALLOCATED(x % start_index)) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 1, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 1
        CALL serialize(io, SIZE(x % start_index, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 1
        CALL serialize(io, LBOUND(x % start_index, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, x % start_index, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    END IF
    CALL serialize(io, '# end_index', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, ALLOCATED(x % end_index), cleanup = .FALSE.)
    IF (ALLOCATED(x % end_index)) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 1, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 1
        CALL serialize(io, SIZE(x % end_index, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 1
        CALL serialize(io, LBOUND(x % end_index, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, x % end_index, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    END IF
    CALL serialize(io, '# start_block', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, ALLOCATED(x % start_block), cleanup = .FALSE.)
    IF (ALLOCATED(x % start_block)) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 1, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 1
        CALL serialize(io, SIZE(x % start_block, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 1
        CALL serialize(io, LBOUND(x % start_block, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, x % start_block, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    END IF
    CALL serialize(io, '# end_block', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, ALLOCATED(x % end_block), cleanup = .FALSE.)
    IF (ALLOCATED(x % end_block)) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 1, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 1
        CALL serialize(io, SIZE(x % end_block, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 1
        CALL serialize(io, LBOUND(x % end_block, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, x % end_block, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    END IF
    IF (nline_local) WRITE(io, '(g0)', ADVANCE = 'no') NEW_LINE('A')
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_t_grid_edges
  SUBROUTINE W_t_grid_vertices(io, x, cleanup, nline)
    USE mo_model_domain, ONLY: t_grid_vertices
    INTEGER :: io
    TYPE(t_grid_vertices), TARGET, INTENT(IN) :: x
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline
    INTEGER :: kmeta, kmeta_0, kmeta_1, kmeta_2, kmeta_3, kmeta_4, kmeta_5, kmeta_6, kmeta_7, kmeta_8, kmeta_9
    LOGICAL :: cleanup_local, nline_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    CALL serialize(io, '# cell_idx', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, ALLOCATED(x % cell_idx), cleanup = .FALSE.)
    IF (ALLOCATED(x % cell_idx)) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 3, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 3
        CALL serialize(io, SIZE(x % cell_idx, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 3
        CALL serialize(io, LBOUND(x % cell_idx, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, x % cell_idx, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    END IF
    CALL serialize(io, '# cell_blk', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, ALLOCATED(x % cell_blk), cleanup = .FALSE.)
    IF (ALLOCATED(x % cell_blk)) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 3, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 3
        CALL serialize(io, SIZE(x % cell_blk, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 3
        CALL serialize(io, LBOUND(x % cell_blk, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, x % cell_blk, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    END IF
    CALL serialize(io, '# edge_idx', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, ALLOCATED(x % edge_idx), cleanup = .FALSE.)
    IF (ALLOCATED(x % edge_idx)) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 3, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 3
        CALL serialize(io, SIZE(x % edge_idx, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 3
        CALL serialize(io, LBOUND(x % edge_idx, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, x % edge_idx, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    END IF
    CALL serialize(io, '# edge_blk', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, ALLOCATED(x % edge_blk), cleanup = .FALSE.)
    IF (ALLOCATED(x % edge_blk)) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 3, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 3
        CALL serialize(io, SIZE(x % edge_blk, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 3
        CALL serialize(io, LBOUND(x % edge_blk, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, x % edge_blk, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    END IF
    CALL serialize(io, '# start_index', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, ALLOCATED(x % start_index), cleanup = .FALSE.)
    IF (ALLOCATED(x % start_index)) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 1, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 1
        CALL serialize(io, SIZE(x % start_index, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 1
        CALL serialize(io, LBOUND(x % start_index, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, x % start_index, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    END IF
    CALL serialize(io, '# end_index', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, ALLOCATED(x % end_index), cleanup = .FALSE.)
    IF (ALLOCATED(x % end_index)) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 1, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 1
        CALL serialize(io, SIZE(x % end_index, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 1
        CALL serialize(io, LBOUND(x % end_index, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, x % end_index, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    END IF
    CALL serialize(io, '# start_block', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, ALLOCATED(x % start_block), cleanup = .FALSE.)
    IF (ALLOCATED(x % start_block)) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 1, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 1
        CALL serialize(io, SIZE(x % start_block, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 1
        CALL serialize(io, LBOUND(x % start_block, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, x % start_block, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    END IF
    CALL serialize(io, '# end_block', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, ALLOCATED(x % end_block), cleanup = .FALSE.)
    IF (ALLOCATED(x % end_block)) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 1, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 1
        CALL serialize(io, SIZE(x % end_block, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 1
        CALL serialize(io, LBOUND(x % end_block, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, x % end_block, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    END IF
    IF (nline_local) WRITE(io, '(g0)', ADVANCE = 'no') NEW_LINE('A')
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_t_grid_vertices
  SUBROUTINE W_t_patch(io, x, cleanup, nline)
    USE mo_model_domain, ONLY: t_patch
    INTEGER :: io
    TYPE(t_patch), TARGET, INTENT(IN) :: x
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline
    INTEGER :: kmeta, kmeta_0, kmeta_1, kmeta_2, kmeta_3, kmeta_4, kmeta_5, kmeta_6, kmeta_7, kmeta_8, kmeta_9
    LOGICAL :: cleanup_local, nline_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    CALL serialize(io, '# nblks_c', cleanup = .FALSE.)
    CALL serialize(io, x % nblks_c, cleanup = .FALSE.)
    CALL serialize(io, '# nblks_e', cleanup = .FALSE.)
    CALL serialize(io, x % nblks_e, cleanup = .FALSE.)
    CALL serialize(io, '# nblks_v', cleanup = .FALSE.)
    CALL serialize(io, x % nblks_v, cleanup = .FALSE.)
    CALL serialize(io, '# cells', cleanup = .FALSE.)
    CALL serialize(io, x % cells, cleanup = .FALSE.)
    CALL serialize(io, '# edges', cleanup = .FALSE.)
    CALL serialize(io, x % edges, cleanup = .FALSE.)
    CALL serialize(io, '# verts', cleanup = .FALSE.)
    CALL serialize(io, x % verts, cleanup = .FALSE.)
    IF (nline_local) WRITE(io, '(g0)', ADVANCE = 'no') NEW_LINE('A')
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_t_patch
  SUBROUTINE W_t_nh_prog(io, x, cleanup, nline)
    USE mo_nonhydro_types, ONLY: t_nh_prog
    INTEGER :: io
    TYPE(t_nh_prog), TARGET, INTENT(IN) :: x
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline
    INTEGER :: kmeta, kmeta_0, kmeta_1, kmeta_2, kmeta_3, kmeta_4, kmeta_5, kmeta_6, kmeta_7, kmeta_8, kmeta_9
    LOGICAL :: cleanup_local, nline_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    CALL serialize(io, '# w', cleanup = .FALSE.)
    CALL serialize(io, '# assoc', cleanup = .FALSE.)
    CALL serialize(io, ASSOCIATED(x % w), cleanup = .FALSE.)
    IF (ASSOCIATED(x % w)) THEN
      kmeta = 0
      CALL serialize(io, "# missing", cleanup = .FALSE.)
      CALL serialize(io, (kmeta == 0), cleanup = .FALSE.)
      CALL serialize(io, x % w, cleanup = .FALSE.)
    END IF
    CALL serialize(io, '# vn', cleanup = .FALSE.)
    CALL serialize(io, '# assoc', cleanup = .FALSE.)
    CALL serialize(io, ASSOCIATED(x % vn), cleanup = .FALSE.)
    IF (ASSOCIATED(x % vn)) THEN
      kmeta = 0
      CALL serialize(io, "# missing", cleanup = .FALSE.)
      CALL serialize(io, (kmeta == 0), cleanup = .FALSE.)
      CALL serialize(io, x % vn, cleanup = .FALSE.)
    END IF
    IF (nline_local) WRITE(io, '(g0)', ADVANCE = 'no') NEW_LINE('A')
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_t_nh_prog
  SUBROUTINE W_t_nh_diag(io, x, cleanup, nline)
    USE mo_nonhydro_types, ONLY: t_nh_diag
    INTEGER :: io
    TYPE(t_nh_diag), TARGET, INTENT(IN) :: x
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline
    INTEGER :: kmeta, kmeta_0, kmeta_1, kmeta_2, kmeta_3, kmeta_4, kmeta_5, kmeta_6, kmeta_7, kmeta_8, kmeta_9
    LOGICAL :: cleanup_local, nline_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    CALL serialize(io, '# vt', cleanup = .FALSE.)
    CALL serialize(io, '# assoc', cleanup = .FALSE.)
    CALL serialize(io, ASSOCIATED(x % vt), cleanup = .FALSE.)
    IF (ASSOCIATED(x % vt)) THEN
      kmeta = 0
      CALL serialize(io, "# missing", cleanup = .FALSE.)
      CALL serialize(io, (kmeta == 0), cleanup = .FALSE.)
      CALL serialize(io, x % vt, cleanup = .FALSE.)
    END IF
    CALL serialize(io, '# vn_ie', cleanup = .FALSE.)
    CALL serialize(io, '# assoc', cleanup = .FALSE.)
    CALL serialize(io, ASSOCIATED(x % vn_ie), cleanup = .FALSE.)
    IF (ASSOCIATED(x % vn_ie)) THEN
      kmeta = 0
      CALL serialize(io, "# missing", cleanup = .FALSE.)
      CALL serialize(io, (kmeta == 0), cleanup = .FALSE.)
      CALL serialize(io, x % vn_ie, cleanup = .FALSE.)
    END IF
    CALL serialize(io, '# w_concorr_c', cleanup = .FALSE.)
    CALL serialize(io, '# assoc', cleanup = .FALSE.)
    CALL serialize(io, ASSOCIATED(x % w_concorr_c), cleanup = .FALSE.)
    IF (ASSOCIATED(x % w_concorr_c)) THEN
      kmeta = 0
      CALL serialize(io, "# missing", cleanup = .FALSE.)
      CALL serialize(io, (kmeta == 0), cleanup = .FALSE.)
      CALL serialize(io, x % w_concorr_c, cleanup = .FALSE.)
    END IF
    CALL serialize(io, '# ddt_vn_apc_pc', cleanup = .FALSE.)
    CALL serialize(io, '# assoc', cleanup = .FALSE.)
    CALL serialize(io, ASSOCIATED(x % ddt_vn_apc_pc), cleanup = .FALSE.)
    IF (ASSOCIATED(x % ddt_vn_apc_pc)) THEN
      kmeta = 0
      CALL serialize(io, "# missing", cleanup = .FALSE.)
      CALL serialize(io, (kmeta == 0), cleanup = .FALSE.)
      CALL serialize(io, x % ddt_vn_apc_pc, cleanup = .FALSE.)
    END IF
    CALL serialize(io, '# ddt_w_adv_pc', cleanup = .FALSE.)
    CALL serialize(io, '# assoc', cleanup = .FALSE.)
    CALL serialize(io, ASSOCIATED(x % ddt_w_adv_pc), cleanup = .FALSE.)
    IF (ASSOCIATED(x % ddt_w_adv_pc)) THEN
      kmeta = 0
      CALL serialize(io, "# missing", cleanup = .FALSE.)
      CALL serialize(io, (kmeta == 0), cleanup = .FALSE.)
      CALL serialize(io, x % ddt_w_adv_pc, cleanup = .FALSE.)
    END IF
    CALL serialize(io, '# max_vcfl_dyn', cleanup = .FALSE.)
    CALL serialize(io, x % max_vcfl_dyn, cleanup = .FALSE.)
    IF (nline_local) WRITE(io, '(g0)', ADVANCE = 'no') NEW_LINE('A')
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_t_nh_diag
  SUBROUTINE W_t_nh_metrics(io, x, cleanup, nline)
    USE mo_nonhydro_types, ONLY: t_nh_metrics
    INTEGER :: io
    TYPE(t_nh_metrics), TARGET, INTENT(IN) :: x
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline
    INTEGER :: kmeta, kmeta_0, kmeta_1, kmeta_2, kmeta_3, kmeta_4, kmeta_5, kmeta_6, kmeta_7, kmeta_8, kmeta_9
    LOGICAL :: cleanup_local, nline_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    CALL serialize(io, '# ddxn_z_full', cleanup = .FALSE.)
    CALL serialize(io, '# assoc', cleanup = .FALSE.)
    CALL serialize(io, ASSOCIATED(x % ddxn_z_full), cleanup = .FALSE.)
    IF (ASSOCIATED(x % ddxn_z_full)) THEN
      kmeta = 0
      CALL serialize(io, "# missing", cleanup = .FALSE.)
      CALL serialize(io, (kmeta == 0), cleanup = .FALSE.)
      CALL serialize(io, x % ddxn_z_full, cleanup = .FALSE.)
    END IF
    CALL serialize(io, '# ddxt_z_full', cleanup = .FALSE.)
    CALL serialize(io, '# assoc', cleanup = .FALSE.)
    CALL serialize(io, ASSOCIATED(x % ddxt_z_full), cleanup = .FALSE.)
    IF (ASSOCIATED(x % ddxt_z_full)) THEN
      kmeta = 0
      CALL serialize(io, "# missing", cleanup = .FALSE.)
      CALL serialize(io, (kmeta == 0), cleanup = .FALSE.)
      CALL serialize(io, x % ddxt_z_full, cleanup = .FALSE.)
    END IF
    CALL serialize(io, '# ddqz_z_full_e', cleanup = .FALSE.)
    CALL serialize(io, '# assoc', cleanup = .FALSE.)
    CALL serialize(io, ASSOCIATED(x % ddqz_z_full_e), cleanup = .FALSE.)
    IF (ASSOCIATED(x % ddqz_z_full_e)) THEN
      kmeta = 0
      CALL serialize(io, "# missing", cleanup = .FALSE.)
      CALL serialize(io, (kmeta == 0), cleanup = .FALSE.)
      CALL serialize(io, x % ddqz_z_full_e, cleanup = .FALSE.)
    END IF
    CALL serialize(io, '# ddqz_z_half', cleanup = .FALSE.)
    CALL serialize(io, '# assoc', cleanup = .FALSE.)
    CALL serialize(io, ASSOCIATED(x % ddqz_z_half), cleanup = .FALSE.)
    IF (ASSOCIATED(x % ddqz_z_half)) THEN
      kmeta = 0
      CALL serialize(io, "# missing", cleanup = .FALSE.)
      CALL serialize(io, (kmeta == 0), cleanup = .FALSE.)
      CALL serialize(io, x % ddqz_z_half, cleanup = .FALSE.)
    END IF
    CALL serialize(io, '# wgtfac_c', cleanup = .FALSE.)
    CALL serialize(io, '# assoc', cleanup = .FALSE.)
    CALL serialize(io, ASSOCIATED(x % wgtfac_c), cleanup = .FALSE.)
    IF (ASSOCIATED(x % wgtfac_c)) THEN
      kmeta = 0
      CALL serialize(io, "# missing", cleanup = .FALSE.)
      CALL serialize(io, (kmeta == 0), cleanup = .FALSE.)
      CALL serialize(io, x % wgtfac_c, cleanup = .FALSE.)
    END IF
    CALL serialize(io, '# wgtfac_e', cleanup = .FALSE.)
    CALL serialize(io, '# assoc', cleanup = .FALSE.)
    CALL serialize(io, ASSOCIATED(x % wgtfac_e), cleanup = .FALSE.)
    IF (ASSOCIATED(x % wgtfac_e)) THEN
      kmeta = 0
      CALL serialize(io, "# missing", cleanup = .FALSE.)
      CALL serialize(io, (kmeta == 0), cleanup = .FALSE.)
      CALL serialize(io, x % wgtfac_e, cleanup = .FALSE.)
    END IF
    CALL serialize(io, '# wgtfacq_e', cleanup = .FALSE.)
    CALL serialize(io, '# assoc', cleanup = .FALSE.)
    CALL serialize(io, ASSOCIATED(x % wgtfacq_e), cleanup = .FALSE.)
    IF (ASSOCIATED(x % wgtfacq_e)) THEN
      kmeta = 0
      CALL serialize(io, "# missing", cleanup = .FALSE.)
      CALL serialize(io, (kmeta == 0), cleanup = .FALSE.)
      CALL serialize(io, x % wgtfacq_e, cleanup = .FALSE.)
    END IF
    CALL serialize(io, '# coeff_gradekin', cleanup = .FALSE.)
    CALL serialize(io, '# assoc', cleanup = .FALSE.)
    CALL serialize(io, ASSOCIATED(x % coeff_gradekin), cleanup = .FALSE.)
    IF (ASSOCIATED(x % coeff_gradekin)) THEN
      kmeta = 0
      CALL serialize(io, "# missing", cleanup = .FALSE.)
      CALL serialize(io, (kmeta == 0), cleanup = .FALSE.)
      CALL serialize(io, x % coeff_gradekin, cleanup = .FALSE.)
    END IF
    CALL serialize(io, '# coeff1_dwdz', cleanup = .FALSE.)
    CALL serialize(io, '# assoc', cleanup = .FALSE.)
    CALL serialize(io, ASSOCIATED(x % coeff1_dwdz), cleanup = .FALSE.)
    IF (ASSOCIATED(x % coeff1_dwdz)) THEN
      kmeta = 0
      CALL serialize(io, "# missing", cleanup = .FALSE.)
      CALL serialize(io, (kmeta == 0), cleanup = .FALSE.)
      CALL serialize(io, x % coeff1_dwdz, cleanup = .FALSE.)
    END IF
    CALL serialize(io, '# coeff2_dwdz', cleanup = .FALSE.)
    CALL serialize(io, '# assoc', cleanup = .FALSE.)
    CALL serialize(io, ASSOCIATED(x % coeff2_dwdz), cleanup = .FALSE.)
    IF (ASSOCIATED(x % coeff2_dwdz)) THEN
      kmeta = 0
      CALL serialize(io, "# missing", cleanup = .FALSE.)
      CALL serialize(io, (kmeta == 0), cleanup = .FALSE.)
      CALL serialize(io, x % coeff2_dwdz, cleanup = .FALSE.)
    END IF
    CALL serialize(io, '# deepatmo_gradh_mc', cleanup = .FALSE.)
    CALL serialize(io, '# assoc', cleanup = .FALSE.)
    CALL serialize(io, ASSOCIATED(x % deepatmo_gradh_mc), cleanup = .FALSE.)
    IF (ASSOCIATED(x % deepatmo_gradh_mc)) THEN
      kmeta = 0
      CALL serialize(io, "# missing", cleanup = .FALSE.)
      CALL serialize(io, (kmeta == 0), cleanup = .FALSE.)
      CALL serialize(io, x % deepatmo_gradh_mc, cleanup = .FALSE.)
    END IF
    CALL serialize(io, '# deepatmo_invr_mc', cleanup = .FALSE.)
    CALL serialize(io, '# assoc', cleanup = .FALSE.)
    CALL serialize(io, ASSOCIATED(x % deepatmo_invr_mc), cleanup = .FALSE.)
    IF (ASSOCIATED(x % deepatmo_invr_mc)) THEN
      kmeta = 0
      CALL serialize(io, "# missing", cleanup = .FALSE.)
      CALL serialize(io, (kmeta == 0), cleanup = .FALSE.)
      CALL serialize(io, x % deepatmo_invr_mc, cleanup = .FALSE.)
    END IF
    CALL serialize(io, '# deepatmo_gradh_ifc', cleanup = .FALSE.)
    CALL serialize(io, '# assoc', cleanup = .FALSE.)
    CALL serialize(io, ASSOCIATED(x % deepatmo_gradh_ifc), cleanup = .FALSE.)
    IF (ASSOCIATED(x % deepatmo_gradh_ifc)) THEN
      kmeta = 0
      CALL serialize(io, "# missing", cleanup = .FALSE.)
      CALL serialize(io, (kmeta == 0), cleanup = .FALSE.)
      CALL serialize(io, x % deepatmo_gradh_ifc, cleanup = .FALSE.)
    END IF
    CALL serialize(io, '# deepatmo_invr_ifc', cleanup = .FALSE.)
    CALL serialize(io, '# assoc', cleanup = .FALSE.)
    CALL serialize(io, ASSOCIATED(x % deepatmo_invr_ifc), cleanup = .FALSE.)
    IF (ASSOCIATED(x % deepatmo_invr_ifc)) THEN
      kmeta = 0
      CALL serialize(io, "# missing", cleanup = .FALSE.)
      CALL serialize(io, (kmeta == 0), cleanup = .FALSE.)
      CALL serialize(io, x % deepatmo_invr_ifc, cleanup = .FALSE.)
    END IF
    IF (nline_local) WRITE(io, '(g0)', ADVANCE = 'no') NEW_LINE('A')
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_t_nh_metrics
  SUBROUTINE W_logical_R_1(io, x, cleanup, nline, meta)
    INTEGER :: io
    LOGICAL, INTENT(IN) :: x(:)
    INTEGER :: k, kmeta, k1
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline, meta
    LOGICAL :: cleanup_local, nline_local, meta_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    meta_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    IF (PRESENT(meta)) meta_local = meta
    IF (meta_local) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 1, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 1
        CALL serialize(io, SIZE(x, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 1
        CALL serialize(io, LBOUND(x, kmeta), cleanup = .FALSE.)
      END DO
    END IF
    CALL serialize(io, "# entries", cleanup = .FALSE.)
    DO k1 = LBOUND(x, 1), UBOUND(x, 1)
      CALL serialize(io, x(k1), cleanup = .FALSE.)
    END DO
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_logical_R_1
  SUBROUTINE W_integer__1_R_1(io, x, cleanup, nline, meta)
    INTEGER :: io
    INTEGER(KIND = 1), INTENT(IN) :: x(:)
    INTEGER :: k, kmeta, k1
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline, meta
    LOGICAL :: cleanup_local, nline_local, meta_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    meta_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    IF (PRESENT(meta)) meta_local = meta
    IF (meta_local) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 1, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 1
        CALL serialize(io, SIZE(x, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 1
        CALL serialize(io, LBOUND(x, kmeta), cleanup = .FALSE.)
      END DO
    END IF
    CALL serialize(io, "# entries", cleanup = .FALSE.)
    DO k1 = LBOUND(x, 1), UBOUND(x, 1)
      CALL serialize(io, x(k1), cleanup = .FALSE.)
    END DO
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_integer__1_R_1
  SUBROUTINE W_integer__2_R_1(io, x, cleanup, nline, meta)
    INTEGER :: io
    INTEGER(KIND = 2), INTENT(IN) :: x(:)
    INTEGER :: k, kmeta, k1
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline, meta
    LOGICAL :: cleanup_local, nline_local, meta_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    meta_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    IF (PRESENT(meta)) meta_local = meta
    IF (meta_local) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 1, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 1
        CALL serialize(io, SIZE(x, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 1
        CALL serialize(io, LBOUND(x, kmeta), cleanup = .FALSE.)
      END DO
    END IF
    CALL serialize(io, "# entries", cleanup = .FALSE.)
    DO k1 = LBOUND(x, 1), UBOUND(x, 1)
      CALL serialize(io, x(k1), cleanup = .FALSE.)
    END DO
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_integer__2_R_1
  SUBROUTINE W_integer__4_R_1(io, x, cleanup, nline, meta)
    INTEGER :: io
    INTEGER(KIND = 4), INTENT(IN) :: x(:)
    INTEGER :: k, kmeta, k1
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline, meta
    LOGICAL :: cleanup_local, nline_local, meta_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    meta_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    IF (PRESENT(meta)) meta_local = meta
    IF (meta_local) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 1, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 1
        CALL serialize(io, SIZE(x, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 1
        CALL serialize(io, LBOUND(x, kmeta), cleanup = .FALSE.)
      END DO
    END IF
    CALL serialize(io, "# entries", cleanup = .FALSE.)
    DO k1 = LBOUND(x, 1), UBOUND(x, 1)
      CALL serialize(io, x(k1), cleanup = .FALSE.)
    END DO
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_integer__4_R_1
  SUBROUTINE W_integer__8_R_1(io, x, cleanup, nline, meta)
    INTEGER :: io
    INTEGER(KIND = 8), INTENT(IN) :: x(:)
    INTEGER :: k, kmeta, k1
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline, meta
    LOGICAL :: cleanup_local, nline_local, meta_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    meta_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    IF (PRESENT(meta)) meta_local = meta
    IF (meta_local) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 1, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 1
        CALL serialize(io, SIZE(x, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 1
        CALL serialize(io, LBOUND(x, kmeta), cleanup = .FALSE.)
      END DO
    END IF
    CALL serialize(io, "# entries", cleanup = .FALSE.)
    DO k1 = LBOUND(x, 1), UBOUND(x, 1)
      CALL serialize(io, x(k1), cleanup = .FALSE.)
    END DO
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_integer__8_R_1
  SUBROUTINE W_real__4_R_1(io, x, cleanup, nline, meta)
    INTEGER :: io
    REAL(KIND = 4), INTENT(IN) :: x(:)
    INTEGER :: k, kmeta, k1
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline, meta
    LOGICAL :: cleanup_local, nline_local, meta_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    meta_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    IF (PRESENT(meta)) meta_local = meta
    IF (meta_local) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 1, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 1
        CALL serialize(io, SIZE(x, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 1
        CALL serialize(io, LBOUND(x, kmeta), cleanup = .FALSE.)
      END DO
    END IF
    CALL serialize(io, "# entries", cleanup = .FALSE.)
    DO k1 = LBOUND(x, 1), UBOUND(x, 1)
      CALL serialize(io, x(k1), cleanup = .FALSE.)
    END DO
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_real__4_R_1
  SUBROUTINE W_real__8_R_1(io, x, cleanup, nline, meta)
    INTEGER :: io
    REAL(KIND = 8), INTENT(IN) :: x(:)
    INTEGER :: k, kmeta, k1
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline, meta
    LOGICAL :: cleanup_local, nline_local, meta_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    meta_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    IF (PRESENT(meta)) meta_local = meta
    IF (meta_local) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 1, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 1
        CALL serialize(io, SIZE(x, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 1
        CALL serialize(io, LBOUND(x, kmeta), cleanup = .FALSE.)
      END DO
    END IF
    CALL serialize(io, "# entries", cleanup = .FALSE.)
    DO k1 = LBOUND(x, 1), UBOUND(x, 1)
      CALL serialize(io, x(k1), cleanup = .FALSE.)
    END DO
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_real__8_R_1
  SUBROUTINE W_logical_R_2(io, x, cleanup, nline, meta)
    INTEGER :: io
    LOGICAL, INTENT(IN) :: x(:, :)
    INTEGER :: k, kmeta, k1, k2
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline, meta
    LOGICAL :: cleanup_local, nline_local, meta_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    meta_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    IF (PRESENT(meta)) meta_local = meta
    IF (meta_local) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 2, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 2
        CALL serialize(io, SIZE(x, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 2
        CALL serialize(io, LBOUND(x, kmeta), cleanup = .FALSE.)
      END DO
    END IF
    CALL serialize(io, "# entries", cleanup = .FALSE.)
    DO k2 = LBOUND(x, 2), UBOUND(x, 2)
      DO k1 = LBOUND(x, 1), UBOUND(x, 1)
        CALL serialize(io, x(k1, k2), cleanup = .FALSE.)
      END DO
    END DO
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_logical_R_2
  SUBROUTINE W_integer__1_R_2(io, x, cleanup, nline, meta)
    INTEGER :: io
    INTEGER(KIND = 1), INTENT(IN) :: x(:, :)
    INTEGER :: k, kmeta, k1, k2
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline, meta
    LOGICAL :: cleanup_local, nline_local, meta_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    meta_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    IF (PRESENT(meta)) meta_local = meta
    IF (meta_local) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 2, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 2
        CALL serialize(io, SIZE(x, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 2
        CALL serialize(io, LBOUND(x, kmeta), cleanup = .FALSE.)
      END DO
    END IF
    CALL serialize(io, "# entries", cleanup = .FALSE.)
    DO k2 = LBOUND(x, 2), UBOUND(x, 2)
      DO k1 = LBOUND(x, 1), UBOUND(x, 1)
        CALL serialize(io, x(k1, k2), cleanup = .FALSE.)
      END DO
    END DO
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_integer__1_R_2
  SUBROUTINE W_integer__2_R_2(io, x, cleanup, nline, meta)
    INTEGER :: io
    INTEGER(KIND = 2), INTENT(IN) :: x(:, :)
    INTEGER :: k, kmeta, k1, k2
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline, meta
    LOGICAL :: cleanup_local, nline_local, meta_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    meta_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    IF (PRESENT(meta)) meta_local = meta
    IF (meta_local) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 2, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 2
        CALL serialize(io, SIZE(x, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 2
        CALL serialize(io, LBOUND(x, kmeta), cleanup = .FALSE.)
      END DO
    END IF
    CALL serialize(io, "# entries", cleanup = .FALSE.)
    DO k2 = LBOUND(x, 2), UBOUND(x, 2)
      DO k1 = LBOUND(x, 1), UBOUND(x, 1)
        CALL serialize(io, x(k1, k2), cleanup = .FALSE.)
      END DO
    END DO
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_integer__2_R_2
  SUBROUTINE W_integer__4_R_2(io, x, cleanup, nline, meta)
    INTEGER :: io
    INTEGER(KIND = 4), INTENT(IN) :: x(:, :)
    INTEGER :: k, kmeta, k1, k2
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline, meta
    LOGICAL :: cleanup_local, nline_local, meta_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    meta_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    IF (PRESENT(meta)) meta_local = meta
    IF (meta_local) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 2, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 2
        CALL serialize(io, SIZE(x, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 2
        CALL serialize(io, LBOUND(x, kmeta), cleanup = .FALSE.)
      END DO
    END IF
    CALL serialize(io, "# entries", cleanup = .FALSE.)
    DO k2 = LBOUND(x, 2), UBOUND(x, 2)
      DO k1 = LBOUND(x, 1), UBOUND(x, 1)
        CALL serialize(io, x(k1, k2), cleanup = .FALSE.)
      END DO
    END DO
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_integer__4_R_2
  SUBROUTINE W_integer__8_R_2(io, x, cleanup, nline, meta)
    INTEGER :: io
    INTEGER(KIND = 8), INTENT(IN) :: x(:, :)
    INTEGER :: k, kmeta, k1, k2
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline, meta
    LOGICAL :: cleanup_local, nline_local, meta_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    meta_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    IF (PRESENT(meta)) meta_local = meta
    IF (meta_local) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 2, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 2
        CALL serialize(io, SIZE(x, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 2
        CALL serialize(io, LBOUND(x, kmeta), cleanup = .FALSE.)
      END DO
    END IF
    CALL serialize(io, "# entries", cleanup = .FALSE.)
    DO k2 = LBOUND(x, 2), UBOUND(x, 2)
      DO k1 = LBOUND(x, 1), UBOUND(x, 1)
        CALL serialize(io, x(k1, k2), cleanup = .FALSE.)
      END DO
    END DO
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_integer__8_R_2
  SUBROUTINE W_real__4_R_2(io, x, cleanup, nline, meta)
    INTEGER :: io
    REAL(KIND = 4), INTENT(IN) :: x(:, :)
    INTEGER :: k, kmeta, k1, k2
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline, meta
    LOGICAL :: cleanup_local, nline_local, meta_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    meta_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    IF (PRESENT(meta)) meta_local = meta
    IF (meta_local) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 2, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 2
        CALL serialize(io, SIZE(x, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 2
        CALL serialize(io, LBOUND(x, kmeta), cleanup = .FALSE.)
      END DO
    END IF
    CALL serialize(io, "# entries", cleanup = .FALSE.)
    DO k2 = LBOUND(x, 2), UBOUND(x, 2)
      DO k1 = LBOUND(x, 1), UBOUND(x, 1)
        CALL serialize(io, x(k1, k2), cleanup = .FALSE.)
      END DO
    END DO
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_real__4_R_2
  SUBROUTINE W_real__8_R_2(io, x, cleanup, nline, meta)
    INTEGER :: io
    REAL(KIND = 8), INTENT(IN) :: x(:, :)
    INTEGER :: k, kmeta, k1, k2
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline, meta
    LOGICAL :: cleanup_local, nline_local, meta_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    meta_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    IF (PRESENT(meta)) meta_local = meta
    IF (meta_local) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 2, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 2
        CALL serialize(io, SIZE(x, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 2
        CALL serialize(io, LBOUND(x, kmeta), cleanup = .FALSE.)
      END DO
    END IF
    CALL serialize(io, "# entries", cleanup = .FALSE.)
    DO k2 = LBOUND(x, 2), UBOUND(x, 2)
      DO k1 = LBOUND(x, 1), UBOUND(x, 1)
        CALL serialize(io, x(k1, k2), cleanup = .FALSE.)
      END DO
    END DO
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_real__8_R_2
  SUBROUTINE W_logical_R_3(io, x, cleanup, nline, meta)
    INTEGER :: io
    LOGICAL, INTENT(IN) :: x(:, :, :)
    INTEGER :: k, kmeta, k1, k2, k3
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline, meta
    LOGICAL :: cleanup_local, nline_local, meta_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    meta_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    IF (PRESENT(meta)) meta_local = meta
    IF (meta_local) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 3, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 3
        CALL serialize(io, SIZE(x, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 3
        CALL serialize(io, LBOUND(x, kmeta), cleanup = .FALSE.)
      END DO
    END IF
    CALL serialize(io, "# entries", cleanup = .FALSE.)
    DO k3 = LBOUND(x, 3), UBOUND(x, 3)
      DO k2 = LBOUND(x, 2), UBOUND(x, 2)
        DO k1 = LBOUND(x, 1), UBOUND(x, 1)
          CALL serialize(io, x(k1, k2, k3), cleanup = .FALSE.)
        END DO
      END DO
    END DO
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_logical_R_3
  SUBROUTINE W_integer__1_R_3(io, x, cleanup, nline, meta)
    INTEGER :: io
    INTEGER(KIND = 1), INTENT(IN) :: x(:, :, :)
    INTEGER :: k, kmeta, k1, k2, k3
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline, meta
    LOGICAL :: cleanup_local, nline_local, meta_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    meta_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    IF (PRESENT(meta)) meta_local = meta
    IF (meta_local) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 3, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 3
        CALL serialize(io, SIZE(x, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 3
        CALL serialize(io, LBOUND(x, kmeta), cleanup = .FALSE.)
      END DO
    END IF
    CALL serialize(io, "# entries", cleanup = .FALSE.)
    DO k3 = LBOUND(x, 3), UBOUND(x, 3)
      DO k2 = LBOUND(x, 2), UBOUND(x, 2)
        DO k1 = LBOUND(x, 1), UBOUND(x, 1)
          CALL serialize(io, x(k1, k2, k3), cleanup = .FALSE.)
        END DO
      END DO
    END DO
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_integer__1_R_3
  SUBROUTINE W_integer__2_R_3(io, x, cleanup, nline, meta)
    INTEGER :: io
    INTEGER(KIND = 2), INTENT(IN) :: x(:, :, :)
    INTEGER :: k, kmeta, k1, k2, k3
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline, meta
    LOGICAL :: cleanup_local, nline_local, meta_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    meta_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    IF (PRESENT(meta)) meta_local = meta
    IF (meta_local) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 3, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 3
        CALL serialize(io, SIZE(x, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 3
        CALL serialize(io, LBOUND(x, kmeta), cleanup = .FALSE.)
      END DO
    END IF
    CALL serialize(io, "# entries", cleanup = .FALSE.)
    DO k3 = LBOUND(x, 3), UBOUND(x, 3)
      DO k2 = LBOUND(x, 2), UBOUND(x, 2)
        DO k1 = LBOUND(x, 1), UBOUND(x, 1)
          CALL serialize(io, x(k1, k2, k3), cleanup = .FALSE.)
        END DO
      END DO
    END DO
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_integer__2_R_3
  SUBROUTINE W_integer__4_R_3(io, x, cleanup, nline, meta)
    INTEGER :: io
    INTEGER(KIND = 4), INTENT(IN) :: x(:, :, :)
    INTEGER :: k, kmeta, k1, k2, k3
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline, meta
    LOGICAL :: cleanup_local, nline_local, meta_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    meta_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    IF (PRESENT(meta)) meta_local = meta
    IF (meta_local) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 3, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 3
        CALL serialize(io, SIZE(x, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 3
        CALL serialize(io, LBOUND(x, kmeta), cleanup = .FALSE.)
      END DO
    END IF
    CALL serialize(io, "# entries", cleanup = .FALSE.)
    DO k3 = LBOUND(x, 3), UBOUND(x, 3)
      DO k2 = LBOUND(x, 2), UBOUND(x, 2)
        DO k1 = LBOUND(x, 1), UBOUND(x, 1)
          CALL serialize(io, x(k1, k2, k3), cleanup = .FALSE.)
        END DO
      END DO
    END DO
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_integer__4_R_3
  SUBROUTINE W_integer__8_R_3(io, x, cleanup, nline, meta)
    INTEGER :: io
    INTEGER(KIND = 8), INTENT(IN) :: x(:, :, :)
    INTEGER :: k, kmeta, k1, k2, k3
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline, meta
    LOGICAL :: cleanup_local, nline_local, meta_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    meta_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    IF (PRESENT(meta)) meta_local = meta
    IF (meta_local) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 3, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 3
        CALL serialize(io, SIZE(x, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 3
        CALL serialize(io, LBOUND(x, kmeta), cleanup = .FALSE.)
      END DO
    END IF
    CALL serialize(io, "# entries", cleanup = .FALSE.)
    DO k3 = LBOUND(x, 3), UBOUND(x, 3)
      DO k2 = LBOUND(x, 2), UBOUND(x, 2)
        DO k1 = LBOUND(x, 1), UBOUND(x, 1)
          CALL serialize(io, x(k1, k2, k3), cleanup = .FALSE.)
        END DO
      END DO
    END DO
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_integer__8_R_3
  SUBROUTINE W_real__4_R_3(io, x, cleanup, nline, meta)
    INTEGER :: io
    REAL(KIND = 4), INTENT(IN) :: x(:, :, :)
    INTEGER :: k, kmeta, k1, k2, k3
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline, meta
    LOGICAL :: cleanup_local, nline_local, meta_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    meta_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    IF (PRESENT(meta)) meta_local = meta
    IF (meta_local) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 3, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 3
        CALL serialize(io, SIZE(x, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 3
        CALL serialize(io, LBOUND(x, kmeta), cleanup = .FALSE.)
      END DO
    END IF
    CALL serialize(io, "# entries", cleanup = .FALSE.)
    DO k3 = LBOUND(x, 3), UBOUND(x, 3)
      DO k2 = LBOUND(x, 2), UBOUND(x, 2)
        DO k1 = LBOUND(x, 1), UBOUND(x, 1)
          CALL serialize(io, x(k1, k2, k3), cleanup = .FALSE.)
        END DO
      END DO
    END DO
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_real__4_R_3
  SUBROUTINE W_real__8_R_3(io, x, cleanup, nline, meta)
    INTEGER :: io
    REAL(KIND = 8), INTENT(IN) :: x(:, :, :)
    INTEGER :: k, kmeta, k1, k2, k3
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline, meta
    LOGICAL :: cleanup_local, nline_local, meta_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    meta_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    IF (PRESENT(meta)) meta_local = meta
    IF (meta_local) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 3, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 3
        CALL serialize(io, SIZE(x, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 3
        CALL serialize(io, LBOUND(x, kmeta), cleanup = .FALSE.)
      END DO
    END IF
    CALL serialize(io, "# entries", cleanup = .FALSE.)
    DO k3 = LBOUND(x, 3), UBOUND(x, 3)
      DO k2 = LBOUND(x, 2), UBOUND(x, 2)
        DO k1 = LBOUND(x, 1), UBOUND(x, 1)
          CALL serialize(io, x(k1, k2, k3), cleanup = .FALSE.)
        END DO
      END DO
    END DO
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_real__8_R_3
  SUBROUTINE W_logical_R_4(io, x, cleanup, nline, meta)
    INTEGER :: io
    LOGICAL, INTENT(IN) :: x(:, :, :, :)
    INTEGER :: k, kmeta, k1, k2, k3, k4
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline, meta
    LOGICAL :: cleanup_local, nline_local, meta_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    meta_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    IF (PRESENT(meta)) meta_local = meta
    IF (meta_local) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 4, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 4
        CALL serialize(io, SIZE(x, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 4
        CALL serialize(io, LBOUND(x, kmeta), cleanup = .FALSE.)
      END DO
    END IF
    CALL serialize(io, "# entries", cleanup = .FALSE.)
    DO k4 = LBOUND(x, 4), UBOUND(x, 4)
      DO k3 = LBOUND(x, 3), UBOUND(x, 3)
        DO k2 = LBOUND(x, 2), UBOUND(x, 2)
          DO k1 = LBOUND(x, 1), UBOUND(x, 1)
            CALL serialize(io, x(k1, k2, k3, k4), cleanup = .FALSE.)
          END DO
        END DO
      END DO
    END DO
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_logical_R_4
  SUBROUTINE W_integer__1_R_4(io, x, cleanup, nline, meta)
    INTEGER :: io
    INTEGER(KIND = 1), INTENT(IN) :: x(:, :, :, :)
    INTEGER :: k, kmeta, k1, k2, k3, k4
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline, meta
    LOGICAL :: cleanup_local, nline_local, meta_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    meta_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    IF (PRESENT(meta)) meta_local = meta
    IF (meta_local) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 4, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 4
        CALL serialize(io, SIZE(x, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 4
        CALL serialize(io, LBOUND(x, kmeta), cleanup = .FALSE.)
      END DO
    END IF
    CALL serialize(io, "# entries", cleanup = .FALSE.)
    DO k4 = LBOUND(x, 4), UBOUND(x, 4)
      DO k3 = LBOUND(x, 3), UBOUND(x, 3)
        DO k2 = LBOUND(x, 2), UBOUND(x, 2)
          DO k1 = LBOUND(x, 1), UBOUND(x, 1)
            CALL serialize(io, x(k1, k2, k3, k4), cleanup = .FALSE.)
          END DO
        END DO
      END DO
    END DO
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_integer__1_R_4
  SUBROUTINE W_integer__2_R_4(io, x, cleanup, nline, meta)
    INTEGER :: io
    INTEGER(KIND = 2), INTENT(IN) :: x(:, :, :, :)
    INTEGER :: k, kmeta, k1, k2, k3, k4
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline, meta
    LOGICAL :: cleanup_local, nline_local, meta_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    meta_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    IF (PRESENT(meta)) meta_local = meta
    IF (meta_local) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 4, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 4
        CALL serialize(io, SIZE(x, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 4
        CALL serialize(io, LBOUND(x, kmeta), cleanup = .FALSE.)
      END DO
    END IF
    CALL serialize(io, "# entries", cleanup = .FALSE.)
    DO k4 = LBOUND(x, 4), UBOUND(x, 4)
      DO k3 = LBOUND(x, 3), UBOUND(x, 3)
        DO k2 = LBOUND(x, 2), UBOUND(x, 2)
          DO k1 = LBOUND(x, 1), UBOUND(x, 1)
            CALL serialize(io, x(k1, k2, k3, k4), cleanup = .FALSE.)
          END DO
        END DO
      END DO
    END DO
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_integer__2_R_4
  SUBROUTINE W_integer__4_R_4(io, x, cleanup, nline, meta)
    INTEGER :: io
    INTEGER(KIND = 4), INTENT(IN) :: x(:, :, :, :)
    INTEGER :: k, kmeta, k1, k2, k3, k4
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline, meta
    LOGICAL :: cleanup_local, nline_local, meta_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    meta_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    IF (PRESENT(meta)) meta_local = meta
    IF (meta_local) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 4, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 4
        CALL serialize(io, SIZE(x, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 4
        CALL serialize(io, LBOUND(x, kmeta), cleanup = .FALSE.)
      END DO
    END IF
    CALL serialize(io, "# entries", cleanup = .FALSE.)
    DO k4 = LBOUND(x, 4), UBOUND(x, 4)
      DO k3 = LBOUND(x, 3), UBOUND(x, 3)
        DO k2 = LBOUND(x, 2), UBOUND(x, 2)
          DO k1 = LBOUND(x, 1), UBOUND(x, 1)
            CALL serialize(io, x(k1, k2, k3, k4), cleanup = .FALSE.)
          END DO
        END DO
      END DO
    END DO
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_integer__4_R_4
  SUBROUTINE W_integer__8_R_4(io, x, cleanup, nline, meta)
    INTEGER :: io
    INTEGER(KIND = 8), INTENT(IN) :: x(:, :, :, :)
    INTEGER :: k, kmeta, k1, k2, k3, k4
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline, meta
    LOGICAL :: cleanup_local, nline_local, meta_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    meta_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    IF (PRESENT(meta)) meta_local = meta
    IF (meta_local) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 4, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 4
        CALL serialize(io, SIZE(x, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 4
        CALL serialize(io, LBOUND(x, kmeta), cleanup = .FALSE.)
      END DO
    END IF
    CALL serialize(io, "# entries", cleanup = .FALSE.)
    DO k4 = LBOUND(x, 4), UBOUND(x, 4)
      DO k3 = LBOUND(x, 3), UBOUND(x, 3)
        DO k2 = LBOUND(x, 2), UBOUND(x, 2)
          DO k1 = LBOUND(x, 1), UBOUND(x, 1)
            CALL serialize(io, x(k1, k2, k3, k4), cleanup = .FALSE.)
          END DO
        END DO
      END DO
    END DO
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_integer__8_R_4
  SUBROUTINE W_real__4_R_4(io, x, cleanup, nline, meta)
    INTEGER :: io
    REAL(KIND = 4), INTENT(IN) :: x(:, :, :, :)
    INTEGER :: k, kmeta, k1, k2, k3, k4
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline, meta
    LOGICAL :: cleanup_local, nline_local, meta_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    meta_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    IF (PRESENT(meta)) meta_local = meta
    IF (meta_local) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 4, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 4
        CALL serialize(io, SIZE(x, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 4
        CALL serialize(io, LBOUND(x, kmeta), cleanup = .FALSE.)
      END DO
    END IF
    CALL serialize(io, "# entries", cleanup = .FALSE.)
    DO k4 = LBOUND(x, 4), UBOUND(x, 4)
      DO k3 = LBOUND(x, 3), UBOUND(x, 3)
        DO k2 = LBOUND(x, 2), UBOUND(x, 2)
          DO k1 = LBOUND(x, 1), UBOUND(x, 1)
            CALL serialize(io, x(k1, k2, k3, k4), cleanup = .FALSE.)
          END DO
        END DO
      END DO
    END DO
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_real__4_R_4
  SUBROUTINE W_real__8_R_4(io, x, cleanup, nline, meta)
    INTEGER :: io
    REAL(KIND = 8), INTENT(IN) :: x(:, :, :, :)
    INTEGER :: k, kmeta, k1, k2, k3, k4
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline, meta
    LOGICAL :: cleanup_local, nline_local, meta_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    meta_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    IF (PRESENT(meta)) meta_local = meta
    IF (meta_local) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 4, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 4
        CALL serialize(io, SIZE(x, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 4
        CALL serialize(io, LBOUND(x, kmeta), cleanup = .FALSE.)
      END DO
    END IF
    CALL serialize(io, "# entries", cleanup = .FALSE.)
    DO k4 = LBOUND(x, 4), UBOUND(x, 4)
      DO k3 = LBOUND(x, 3), UBOUND(x, 3)
        DO k2 = LBOUND(x, 2), UBOUND(x, 2)
          DO k1 = LBOUND(x, 1), UBOUND(x, 1)
            CALL serialize(io, x(k1, k2, k3, k4), cleanup = .FALSE.)
          END DO
        END DO
      END DO
    END DO
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_real__8_R_4
  SUBROUTINE W_logical(io, x, cleanup, nline)
    CHARACTER(LEN = 50) :: buf
    INTEGER :: io
    LOGICAL, INTENT(IN) :: x
    INTEGER :: y
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline
    LOGICAL :: cleanup_local, nline_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    y = MERGE(1, 0, x)
    WRITE(io, '(g0)', ADVANCE = 'no') y
    IF (nline_local) WRITE(io, '(g0)', ADVANCE = 'no') NEW_LINE('A')
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_logical
  SUBROUTINE W_integer1(io, x, cleanup, nline)
    CHARACTER(LEN = 50) :: buf
    INTEGER :: io
    INTEGER(KIND = 1), INTENT(IN) :: x
    INTEGER :: y
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline
    LOGICAL :: cleanup_local, nline_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    WRITE(io, '(g0)', ADVANCE = 'no') x
    IF (nline_local) WRITE(io, '(g0)', ADVANCE = 'no') NEW_LINE('A')
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_integer1
  SUBROUTINE W_integer2(io, x, cleanup, nline)
    CHARACTER(LEN = 50) :: buf
    INTEGER :: io
    INTEGER(KIND = 2), INTENT(IN) :: x
    INTEGER :: y
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline
    LOGICAL :: cleanup_local, nline_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    WRITE(io, '(g0)', ADVANCE = 'no') x
    IF (nline_local) WRITE(io, '(g0)', ADVANCE = 'no') NEW_LINE('A')
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_integer2
  SUBROUTINE W_integer4(io, x, cleanup, nline)
    CHARACTER(LEN = 50) :: buf
    INTEGER :: io
    INTEGER(KIND = 4), INTENT(IN) :: x
    INTEGER :: y
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline
    LOGICAL :: cleanup_local, nline_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    WRITE(io, '(g0)', ADVANCE = 'no') x
    IF (nline_local) WRITE(io, '(g0)', ADVANCE = 'no') NEW_LINE('A')
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_integer4
  SUBROUTINE W_integer8(io, x, cleanup, nline)
    CHARACTER(LEN = 50) :: buf
    INTEGER :: io
    INTEGER(KIND = 8), INTENT(IN) :: x
    INTEGER :: y
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline
    LOGICAL :: cleanup_local, nline_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    WRITE(io, '(g0)', ADVANCE = 'no') x
    IF (nline_local) WRITE(io, '(g0)', ADVANCE = 'no') NEW_LINE('A')
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_integer8
  SUBROUTINE W_real4(io, x, cleanup, nline)
    CHARACTER(LEN = 50) :: buf
    INTEGER :: io
    REAL(KIND = 4), INTENT(IN) :: x
    INTEGER :: y
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline
    LOGICAL :: cleanup_local, nline_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    WRITE(buf, '(e28.20)') x
    WRITE(io, '(A)', ADVANCE = 'no') TRIM(ADJUSTL(buf))
    IF (nline_local) WRITE(io, '(g0)', ADVANCE = 'no') NEW_LINE('A')
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_real4
  SUBROUTINE W_real8(io, x, cleanup, nline)
    CHARACTER(LEN = 50) :: buf
    INTEGER :: io
    REAL(KIND = 8), INTENT(IN) :: x
    INTEGER :: y
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline
    LOGICAL :: cleanup_local, nline_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    WRITE(buf, '(e28.20)') x
    WRITE(io, '(A)', ADVANCE = 'no') TRIM(ADJUSTL(buf))
    IF (nline_local) WRITE(io, '(g0)', ADVANCE = 'no') NEW_LINE('A')
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_real8
  SUBROUTINE serialize_global_data(io)
    USE mo_init_vgrid, ONLY: nflatlev => nflatlev
    USE mo_mpi, ONLY: i_am_accel_node => i_am_accel_node
    USE mo_nonhydrostatic_config, ONLY: lextra_diffu => lextra_diffu
    USE mo_parallel_config, ONLY: nproma => nproma
    USE mo_run_config, ONLY: timers_level => timers_level
    USE mo_timer, ONLY: timer_solve_nh_veltend => timer_solve_nh_veltend
    USE mo_timer, ONLY: timer_intp => timer_intp
    USE mo_vertical_grid, ONLY: nrdmax => nrdmax
    INTEGER :: io
    CALL serialize(io, "# nflatlev", cleanup = .FALSE.)
    CALL serialize(io, nflatlev, cleanup = .FALSE.)
    CALL serialize(io, "# i_am_accel_node", cleanup = .FALSE.)
    CALL serialize(io, i_am_accel_node, cleanup = .FALSE.)
    CALL serialize(io, "# lextra_diffu", cleanup = .FALSE.)
    CALL serialize(io, lextra_diffu, cleanup = .FALSE.)
    CALL serialize(io, "# nproma", cleanup = .FALSE.)
    CALL serialize(io, nproma, cleanup = .FALSE.)
    CALL serialize(io, "# timers_level", cleanup = .FALSE.)
    CALL serialize(io, timers_level, cleanup = .FALSE.)
    CALL serialize(io, "# timer_solve_nh_veltend", cleanup = .FALSE.)
    CALL serialize(io, timer_solve_nh_veltend, cleanup = .FALSE.)
    CALL serialize(io, "# timer_intp", cleanup = .FALSE.)
    CALL serialize(io, timer_intp, cleanup = .FALSE.)
    CALL serialize(io, "# nrdmax", cleanup = .FALSE.)
    CALL serialize(io, nrdmax, cleanup = .FALSE.)
    CLOSE(UNIT = io)
  END SUBROUTINE serialize_global_data
END MODULE serde