MODULE serdesl
  IMPLICIT NONE
  INTERFACE serialize
    MODULE PROCEDURE :: character_2s
    MODULE PROCEDURE cloud_type_2s, pdf_sampler_type_2s, config_type_2s, flux_type_2s, randomnumberstream_2s, single_level_type_2s, logical_2s1, integer__1_2s1, integer__2_2s1, integer__4_2s1, integer__8_2s1, real__4_2s1, real__8_2s1, logical_2s2, integer__1_2s2, integer__2_2s2, integer__4_2s2, integer__8_2s2, real__4_2s2, real__8_2s2, logical_2s3, integer__1_2s3, integer__2_2s3, integer__4_2s3, integer__8_2s3, real__4_2s3, real__8_2s3, logical_2s4, integer__1_2s4, integer__2_2s4, integer__4_2s4, integer__8_2s4, real__4_2s4, real__8_2s4, logical_2s, integer1_2s, integer2_2s, integer4_2s, integer8_2s, real4_2s, real8_2s
  END INTERFACE serialize
  CONTAINS
  SUBROUTINE write_to(path, s)
    CHARACTER(LEN = *), INTENT(IN) :: path
    CHARACTER(LEN = *), INTENT(IN) :: s
    INTEGER :: io
    OPEN(NEWUNIT = io, FILE = path, STATUS = "replace", ACTION = "write")
    WRITE(io, *) s
    CLOSE(UNIT = io)
  END SUBROUTINE write_to
  FUNCTION add_line(r, l) RESULT(s)
    CHARACTER(LEN = *), INTENT(IN) :: r
    CHARACTER(LEN = *), INTENT(IN) :: l
    CHARACTER(LEN = :), ALLOCATABLE :: s
    s = r // TRIM(l) // NEW_LINE('A')
  END FUNCTION add_line
  FUNCTION character_2s(x) RESULT(s)
    CHARACTER(LEN = *), INTENT(IN) :: x
    CHARACTER(LEN = :), ALLOCATABLE :: s
    ALLOCATE(CHARACTER(LEN = LEN(x) + 1)::s)
    WRITE(s, '(g0)') TRIM(x)
    s = TRIM(s)
  END FUNCTION character_2s
  FUNCTION cloud_type_2s(x) RESULT(s)
    USE radiation_cloud, ONLY: cloud_type
    TYPE(cloud_type), TARGET, INTENT(IN) :: x
    CHARACTER(LEN = :), ALLOCATABLE :: s
    INTEGER :: kmeta, kmeta_0, kmeta_1, kmeta_2, kmeta_3, kmeta_4, kmeta_5, kmeta_6, kmeta_7, kmeta_8, kmeta_9
    s = ""
    s = add_line(s, '# fraction')
    s = add_line(s, '# alloc')
    s = add_line(s, serialize(ALLOCATED(x % fraction)))
    IF (ALLOCATED(x % fraction)) THEN
      s = add_line(s, "# rank")
      s = add_line(s, serialize(2))
      s = add_line(s, "# size")
      DO kmeta = 1, 2
        s = add_line(s, serialize(SIZE(x % fraction, kmeta)))
      END DO
      s = add_line(s, "# lbound")
      DO kmeta = 1, 2
        s = add_line(s, serialize(LBOUND(x % fraction, kmeta)))
      END DO
      s = add_line(s, serialize(x % fraction))
    END IF
    s = add_line(s, '# fractional_std')
    s = add_line(s, '# alloc')
    s = add_line(s, serialize(ALLOCATED(x % fractional_std)))
    IF (ALLOCATED(x % fractional_std)) THEN
      s = add_line(s, "# rank")
      s = add_line(s, serialize(2))
      s = add_line(s, "# size")
      DO kmeta = 1, 2
        s = add_line(s, serialize(SIZE(x % fractional_std, kmeta)))
      END DO
      s = add_line(s, "# lbound")
      DO kmeta = 1, 2
        s = add_line(s, serialize(LBOUND(x % fractional_std, kmeta)))
      END DO
      s = add_line(s, serialize(x % fractional_std))
    END IF
    s = add_line(s, '# overlap_param')
    s = add_line(s, '# alloc')
    s = add_line(s, serialize(ALLOCATED(x % overlap_param)))
    IF (ALLOCATED(x % overlap_param)) THEN
      s = add_line(s, "# rank")
      s = add_line(s, serialize(2))
      s = add_line(s, "# size")
      DO kmeta = 1, 2
        s = add_line(s, serialize(SIZE(x % overlap_param, kmeta)))
      END DO
      s = add_line(s, "# lbound")
      DO kmeta = 1, 2
        s = add_line(s, serialize(LBOUND(x % overlap_param, kmeta)))
      END DO
      s = add_line(s, serialize(x % overlap_param))
    END IF
    IF (LEN(s) > 0) s = s(: LEN(s) - 1)
  END FUNCTION cloud_type_2s
  FUNCTION pdf_sampler_type_2s(x) RESULT(s)
    USE radiation_pdf_sampler, ONLY: pdf_sampler_type
    TYPE(pdf_sampler_type), TARGET, INTENT(IN) :: x
    CHARACTER(LEN = :), ALLOCATABLE :: s
    INTEGER :: kmeta, kmeta_0, kmeta_1, kmeta_2, kmeta_3, kmeta_4, kmeta_5, kmeta_6, kmeta_7, kmeta_8, kmeta_9
    s = ""
    s = add_line(s, '# ncdf')
    s = add_line(s, serialize(x % ncdf))
    s = add_line(s, '# nfsd')
    s = add_line(s, serialize(x % nfsd))
    s = add_line(s, '# fsd1')
    s = add_line(s, serialize(x % fsd1))
    s = add_line(s, '# inv_fsd_interval')
    s = add_line(s, serialize(x % inv_fsd_interval))
    s = add_line(s, '# val')
    s = add_line(s, '# alloc')
    s = add_line(s, serialize(ALLOCATED(x % val)))
    IF (ALLOCATED(x % val)) THEN
      s = add_line(s, "# rank")
      s = add_line(s, serialize(2))
      s = add_line(s, "# size")
      DO kmeta = 1, 2
        s = add_line(s, serialize(SIZE(x % val, kmeta)))
      END DO
      s = add_line(s, "# lbound")
      DO kmeta = 1, 2
        s = add_line(s, serialize(LBOUND(x % val, kmeta)))
      END DO
      s = add_line(s, serialize(x % val))
    END IF
    IF (LEN(s) > 0) s = s(: LEN(s) - 1)
  END FUNCTION pdf_sampler_type_2s
  FUNCTION config_type_2s(x) RESULT(s)
    USE radiation_config, ONLY: config_type
    TYPE(config_type), TARGET, INTENT(IN) :: x
    CHARACTER(LEN = :), ALLOCATABLE :: s
    INTEGER :: kmeta, kmeta_0, kmeta_1, kmeta_2, kmeta_3, kmeta_4, kmeta_5, kmeta_6, kmeta_7, kmeta_8, kmeta_9
    s = ""
    s = add_line(s, '# i_band_from_reordered_g_lw')
    s = add_line(s, '# alloc')
    s = add_line(s, serialize(ALLOCATED(x % i_band_from_reordered_g_lw)))
    IF (ALLOCATED(x % i_band_from_reordered_g_lw)) THEN
      s = add_line(s, "# rank")
      s = add_line(s, serialize(1))
      s = add_line(s, "# size")
      DO kmeta = 1, 1
        s = add_line(s, serialize(SIZE(x % i_band_from_reordered_g_lw, kmeta)))
      END DO
      s = add_line(s, "# lbound")
      DO kmeta = 1, 1
        s = add_line(s, serialize(LBOUND(x % i_band_from_reordered_g_lw, kmeta)))
      END DO
      s = add_line(s, serialize(x % i_band_from_reordered_g_lw))
    END IF
    s = add_line(s, '# pdf_sampler')
    s = add_line(s, serialize(x % pdf_sampler))
    IF (LEN(s) > 0) s = s(: LEN(s) - 1)
  END FUNCTION config_type_2s
  FUNCTION flux_type_2s(x) RESULT(s)
    USE radiation_flux, ONLY: flux_type
    TYPE(flux_type), TARGET, INTENT(IN) :: x
    CHARACTER(LEN = :), ALLOCATABLE :: s
    INTEGER :: kmeta, kmeta_0, kmeta_1, kmeta_2, kmeta_3, kmeta_4, kmeta_5, kmeta_6, kmeta_7, kmeta_8, kmeta_9
    s = ""
    s = add_line(s, '# lw_up')
    s = add_line(s, '# alloc')
    s = add_line(s, serialize(ALLOCATED(x % lw_up)))
    IF (ALLOCATED(x % lw_up)) THEN
      s = add_line(s, "# rank")
      s = add_line(s, serialize(2))
      s = add_line(s, "# size")
      DO kmeta = 1, 2
        s = add_line(s, serialize(SIZE(x % lw_up, kmeta)))
      END DO
      s = add_line(s, "# lbound")
      DO kmeta = 1, 2
        s = add_line(s, serialize(LBOUND(x % lw_up, kmeta)))
      END DO
      s = add_line(s, serialize(x % lw_up))
    END IF
    s = add_line(s, '# lw_dn')
    s = add_line(s, '# alloc')
    s = add_line(s, serialize(ALLOCATED(x % lw_dn)))
    IF (ALLOCATED(x % lw_dn)) THEN
      s = add_line(s, "# rank")
      s = add_line(s, serialize(2))
      s = add_line(s, "# size")
      DO kmeta = 1, 2
        s = add_line(s, serialize(SIZE(x % lw_dn, kmeta)))
      END DO
      s = add_line(s, "# lbound")
      DO kmeta = 1, 2
        s = add_line(s, serialize(LBOUND(x % lw_dn, kmeta)))
      END DO
      s = add_line(s, serialize(x % lw_dn))
    END IF
    s = add_line(s, '# lw_up_clear')
    s = add_line(s, '# alloc')
    s = add_line(s, serialize(ALLOCATED(x % lw_up_clear)))
    IF (ALLOCATED(x % lw_up_clear)) THEN
      s = add_line(s, "# rank")
      s = add_line(s, serialize(2))
      s = add_line(s, "# size")
      DO kmeta = 1, 2
        s = add_line(s, serialize(SIZE(x % lw_up_clear, kmeta)))
      END DO
      s = add_line(s, "# lbound")
      DO kmeta = 1, 2
        s = add_line(s, serialize(LBOUND(x % lw_up_clear, kmeta)))
      END DO
      s = add_line(s, serialize(x % lw_up_clear))
    END IF
    s = add_line(s, '# lw_dn_clear')
    s = add_line(s, '# alloc')
    s = add_line(s, serialize(ALLOCATED(x % lw_dn_clear)))
    IF (ALLOCATED(x % lw_dn_clear)) THEN
      s = add_line(s, "# rank")
      s = add_line(s, serialize(2))
      s = add_line(s, "# size")
      DO kmeta = 1, 2
        s = add_line(s, serialize(SIZE(x % lw_dn_clear, kmeta)))
      END DO
      s = add_line(s, "# lbound")
      DO kmeta = 1, 2
        s = add_line(s, serialize(LBOUND(x % lw_dn_clear, kmeta)))
      END DO
      s = add_line(s, serialize(x % lw_dn_clear))
    END IF
    s = add_line(s, '# lw_dn_surf_g')
    s = add_line(s, '# alloc')
    s = add_line(s, serialize(ALLOCATED(x % lw_dn_surf_g)))
    IF (ALLOCATED(x % lw_dn_surf_g)) THEN
      s = add_line(s, "# rank")
      s = add_line(s, serialize(2))
      s = add_line(s, "# size")
      DO kmeta = 1, 2
        s = add_line(s, serialize(SIZE(x % lw_dn_surf_g, kmeta)))
      END DO
      s = add_line(s, "# lbound")
      DO kmeta = 1, 2
        s = add_line(s, serialize(LBOUND(x % lw_dn_surf_g, kmeta)))
      END DO
      s = add_line(s, serialize(x % lw_dn_surf_g))
    END IF
    s = add_line(s, '# lw_dn_surf_clear_g')
    s = add_line(s, '# alloc')
    s = add_line(s, serialize(ALLOCATED(x % lw_dn_surf_clear_g)))
    IF (ALLOCATED(x % lw_dn_surf_clear_g)) THEN
      s = add_line(s, "# rank")
      s = add_line(s, serialize(2))
      s = add_line(s, "# size")
      DO kmeta = 1, 2
        s = add_line(s, serialize(SIZE(x % lw_dn_surf_clear_g, kmeta)))
      END DO
      s = add_line(s, "# lbound")
      DO kmeta = 1, 2
        s = add_line(s, serialize(LBOUND(x % lw_dn_surf_clear_g, kmeta)))
      END DO
      s = add_line(s, serialize(x % lw_dn_surf_clear_g))
    END IF
    s = add_line(s, '# cloud_cover_lw')
    s = add_line(s, '# alloc')
    s = add_line(s, serialize(ALLOCATED(x % cloud_cover_lw)))
    IF (ALLOCATED(x % cloud_cover_lw)) THEN
      s = add_line(s, "# rank")
      s = add_line(s, serialize(1))
      s = add_line(s, "# size")
      DO kmeta = 1, 1
        s = add_line(s, serialize(SIZE(x % cloud_cover_lw, kmeta)))
      END DO
      s = add_line(s, "# lbound")
      DO kmeta = 1, 1
        s = add_line(s, serialize(LBOUND(x % cloud_cover_lw, kmeta)))
      END DO
      s = add_line(s, serialize(x % cloud_cover_lw))
    END IF
    IF (LEN(s) > 0) s = s(: LEN(s) - 1)
  END FUNCTION flux_type_2s
  FUNCTION randomnumberstream_2s(x) RESULT(s)
    USE random_numbers_mix, ONLY: randomnumberstream
    TYPE(randomnumberstream), TARGET, INTENT(IN) :: x
    CHARACTER(LEN = :), ALLOCATABLE :: s
    INTEGER :: kmeta, kmeta_0, kmeta_1, kmeta_2, kmeta_3, kmeta_4, kmeta_5, kmeta_6, kmeta_7, kmeta_8, kmeta_9
    s = ""
    s = add_line(s, '# iused')
    s = add_line(s, serialize(x % iused))
    s = add_line(s, '# inittest')
    s = add_line(s, serialize(x % inittest))
    s = add_line(s, '# ix')
    s = add_line(s, "# rank")
    s = add_line(s, serialize(1))
    s = add_line(s, "# size")
    DO kmeta = 1, 1
      s = add_line(s, serialize(SIZE(x % ix, kmeta)))
    END DO
    s = add_line(s, "# lbound")
    DO kmeta = 1, 1
      s = add_line(s, serialize(LBOUND(x % ix, kmeta)))
    END DO
    s = add_line(s, serialize(x % ix))
    s = add_line(s, '# zrm')
    s = add_line(s, serialize(x % zrm))
    IF (LEN(s) > 0) s = s(: LEN(s) - 1)
  END FUNCTION randomnumberstream_2s
  FUNCTION single_level_type_2s(x) RESULT(s)
    USE radiation_single_level, ONLY: single_level_type
    TYPE(single_level_type), TARGET, INTENT(IN) :: x
    CHARACTER(LEN = :), ALLOCATABLE :: s
    INTEGER :: kmeta, kmeta_0, kmeta_1, kmeta_2, kmeta_3, kmeta_4, kmeta_5, kmeta_6, kmeta_7, kmeta_8, kmeta_9
    s = ""
    s = add_line(s, '# iseed')
    s = add_line(s, '# alloc')
    s = add_line(s, serialize(ALLOCATED(x % iseed)))
    IF (ALLOCATED(x % iseed)) THEN
      s = add_line(s, "# rank")
      s = add_line(s, serialize(1))
      s = add_line(s, "# size")
      DO kmeta = 1, 1
        s = add_line(s, serialize(SIZE(x % iseed, kmeta)))
      END DO
      s = add_line(s, "# lbound")
      DO kmeta = 1, 1
        s = add_line(s, serialize(LBOUND(x % iseed, kmeta)))
      END DO
      s = add_line(s, serialize(x % iseed))
    END IF
    IF (LEN(s) > 0) s = s(: LEN(s) - 1)
  END FUNCTION single_level_type_2s
  FUNCTION logical_2s1(a) RESULT(s)
    LOGICAL, INTENT(IN) :: a(:)
    CHARACTER(LEN = :), ALLOCATABLE :: s
    INTEGER :: k, k1
    s = ""
    s = add_line(s, "# entries")
    DO k1 = LBOUND(a, 1), UBOUND(a, 1)
      s = add_line(s, serialize(a(k1)))
    END DO
    IF (LEN(s) > 0) s = s(: LEN(s) - 1)
  END FUNCTION logical_2s1
  FUNCTION integer__1_2s1(a) RESULT(s)
    INTEGER(KIND = 1), INTENT(IN) :: a(:)
    CHARACTER(LEN = :), ALLOCATABLE :: s
    INTEGER :: k, k1
    s = ""
    s = add_line(s, "# entries")
    DO k1 = LBOUND(a, 1), UBOUND(a, 1)
      s = add_line(s, serialize(a(k1)))
    END DO
    IF (LEN(s) > 0) s = s(: LEN(s) - 1)
  END FUNCTION integer__1_2s1
  FUNCTION integer__2_2s1(a) RESULT(s)
    INTEGER(KIND = 2), INTENT(IN) :: a(:)
    CHARACTER(LEN = :), ALLOCATABLE :: s
    INTEGER :: k, k1
    s = ""
    s = add_line(s, "# entries")
    DO k1 = LBOUND(a, 1), UBOUND(a, 1)
      s = add_line(s, serialize(a(k1)))
    END DO
    IF (LEN(s) > 0) s = s(: LEN(s) - 1)
  END FUNCTION integer__2_2s1
  FUNCTION integer__4_2s1(a) RESULT(s)
    INTEGER(KIND = 4), INTENT(IN) :: a(:)
    CHARACTER(LEN = :), ALLOCATABLE :: s
    INTEGER :: k, k1
    s = ""
    s = add_line(s, "# entries")
    DO k1 = LBOUND(a, 1), UBOUND(a, 1)
      s = add_line(s, serialize(a(k1)))
    END DO
    IF (LEN(s) > 0) s = s(: LEN(s) - 1)
  END FUNCTION integer__4_2s1
  FUNCTION integer__8_2s1(a) RESULT(s)
    INTEGER(KIND = 8), INTENT(IN) :: a(:)
    CHARACTER(LEN = :), ALLOCATABLE :: s
    INTEGER :: k, k1
    s = ""
    s = add_line(s, "# entries")
    DO k1 = LBOUND(a, 1), UBOUND(a, 1)
      s = add_line(s, serialize(a(k1)))
    END DO
    IF (LEN(s) > 0) s = s(: LEN(s) - 1)
  END FUNCTION integer__8_2s1
  FUNCTION real__4_2s1(a) RESULT(s)
    REAL(KIND = 4), INTENT(IN) :: a(:)
    CHARACTER(LEN = :), ALLOCATABLE :: s
    INTEGER :: k, k1
    s = ""
    s = add_line(s, "# entries")
    DO k1 = LBOUND(a, 1), UBOUND(a, 1)
      s = add_line(s, serialize(a(k1)))
    END DO
    IF (LEN(s) > 0) s = s(: LEN(s) - 1)
  END FUNCTION real__4_2s1
  FUNCTION real__8_2s1(a) RESULT(s)
    REAL(KIND = 8), INTENT(IN) :: a(:)
    CHARACTER(LEN = :), ALLOCATABLE :: s
    INTEGER :: k, k1
    s = ""
    s = add_line(s, "# entries")
    DO k1 = LBOUND(a, 1), UBOUND(a, 1)
      s = add_line(s, serialize(a(k1)))
    END DO
    IF (LEN(s) > 0) s = s(: LEN(s) - 1)
  END FUNCTION real__8_2s1
  FUNCTION logical_2s2(a) RESULT(s)
    LOGICAL, INTENT(IN) :: a(:, :)
    CHARACTER(LEN = :), ALLOCATABLE :: s
    INTEGER :: k, k1, k2
    s = ""
    s = add_line(s, "# entries")
    DO k1 = LBOUND(a, 1), UBOUND(a, 1)
      DO k2 = LBOUND(a, 2), UBOUND(a, 2)
        s = add_line(s, serialize(a(k1, k2)))
      END DO
    END DO
    IF (LEN(s) > 0) s = s(: LEN(s) - 1)
  END FUNCTION logical_2s2
  FUNCTION integer__1_2s2(a) RESULT(s)
    INTEGER(KIND = 1), INTENT(IN) :: a(:, :)
    CHARACTER(LEN = :), ALLOCATABLE :: s
    INTEGER :: k, k1, k2
    s = ""
    s = add_line(s, "# entries")
    DO k1 = LBOUND(a, 1), UBOUND(a, 1)
      DO k2 = LBOUND(a, 2), UBOUND(a, 2)
        s = add_line(s, serialize(a(k1, k2)))
      END DO
    END DO
    IF (LEN(s) > 0) s = s(: LEN(s) - 1)
  END FUNCTION integer__1_2s2
  FUNCTION integer__2_2s2(a) RESULT(s)
    INTEGER(KIND = 2), INTENT(IN) :: a(:, :)
    CHARACTER(LEN = :), ALLOCATABLE :: s
    INTEGER :: k, k1, k2
    s = ""
    s = add_line(s, "# entries")
    DO k1 = LBOUND(a, 1), UBOUND(a, 1)
      DO k2 = LBOUND(a, 2), UBOUND(a, 2)
        s = add_line(s, serialize(a(k1, k2)))
      END DO
    END DO
    IF (LEN(s) > 0) s = s(: LEN(s) - 1)
  END FUNCTION integer__2_2s2
  FUNCTION integer__4_2s2(a) RESULT(s)
    INTEGER(KIND = 4), INTENT(IN) :: a(:, :)
    CHARACTER(LEN = :), ALLOCATABLE :: s
    INTEGER :: k, k1, k2
    s = ""
    s = add_line(s, "# entries")
    DO k1 = LBOUND(a, 1), UBOUND(a, 1)
      DO k2 = LBOUND(a, 2), UBOUND(a, 2)
        s = add_line(s, serialize(a(k1, k2)))
      END DO
    END DO
    IF (LEN(s) > 0) s = s(: LEN(s) - 1)
  END FUNCTION integer__4_2s2
  FUNCTION integer__8_2s2(a) RESULT(s)
    INTEGER(KIND = 8), INTENT(IN) :: a(:, :)
    CHARACTER(LEN = :), ALLOCATABLE :: s
    INTEGER :: k, k1, k2
    s = ""
    s = add_line(s, "# entries")
    DO k1 = LBOUND(a, 1), UBOUND(a, 1)
      DO k2 = LBOUND(a, 2), UBOUND(a, 2)
        s = add_line(s, serialize(a(k1, k2)))
      END DO
    END DO
    IF (LEN(s) > 0) s = s(: LEN(s) - 1)
  END FUNCTION integer__8_2s2
  FUNCTION real__4_2s2(a) RESULT(s)
    REAL(KIND = 4), INTENT(IN) :: a(:, :)
    CHARACTER(LEN = :), ALLOCATABLE :: s
    INTEGER :: k, k1, k2
    s = ""
    s = add_line(s, "# entries")
    DO k1 = LBOUND(a, 1), UBOUND(a, 1)
      DO k2 = LBOUND(a, 2), UBOUND(a, 2)
        s = add_line(s, serialize(a(k1, k2)))
      END DO
    END DO
    IF (LEN(s) > 0) s = s(: LEN(s) - 1)
  END FUNCTION real__4_2s2
  FUNCTION real__8_2s2(a) RESULT(s)
    REAL(KIND = 8), INTENT(IN) :: a(:, :)
    CHARACTER(LEN = :), ALLOCATABLE :: s
    INTEGER :: k, k1, k2
    s = ""
    s = add_line(s, "# entries")
    DO k1 = LBOUND(a, 1), UBOUND(a, 1)
      DO k2 = LBOUND(a, 2), UBOUND(a, 2)
        s = add_line(s, serialize(a(k1, k2)))
      END DO
    END DO
    IF (LEN(s) > 0) s = s(: LEN(s) - 1)
  END FUNCTION real__8_2s2
  FUNCTION logical_2s3(a) RESULT(s)
    LOGICAL, INTENT(IN) :: a(:, :, :)
    CHARACTER(LEN = :), ALLOCATABLE :: s
    INTEGER :: k, k1, k2, k3
    s = ""
    s = add_line(s, "# entries")
    DO k1 = LBOUND(a, 1), UBOUND(a, 1)
      DO k2 = LBOUND(a, 2), UBOUND(a, 2)
        DO k3 = LBOUND(a, 3), UBOUND(a, 3)
          s = add_line(s, serialize(a(k1, k2, k3)))
        END DO
      END DO
    END DO
    IF (LEN(s) > 0) s = s(: LEN(s) - 1)
  END FUNCTION logical_2s3
  FUNCTION integer__1_2s3(a) RESULT(s)
    INTEGER(KIND = 1), INTENT(IN) :: a(:, :, :)
    CHARACTER(LEN = :), ALLOCATABLE :: s
    INTEGER :: k, k1, k2, k3
    s = ""
    s = add_line(s, "# entries")
    DO k1 = LBOUND(a, 1), UBOUND(a, 1)
      DO k2 = LBOUND(a, 2), UBOUND(a, 2)
        DO k3 = LBOUND(a, 3), UBOUND(a, 3)
          s = add_line(s, serialize(a(k1, k2, k3)))
        END DO
      END DO
    END DO
    IF (LEN(s) > 0) s = s(: LEN(s) - 1)
  END FUNCTION integer__1_2s3
  FUNCTION integer__2_2s3(a) RESULT(s)
    INTEGER(KIND = 2), INTENT(IN) :: a(:, :, :)
    CHARACTER(LEN = :), ALLOCATABLE :: s
    INTEGER :: k, k1, k2, k3
    s = ""
    s = add_line(s, "# entries")
    DO k1 = LBOUND(a, 1), UBOUND(a, 1)
      DO k2 = LBOUND(a, 2), UBOUND(a, 2)
        DO k3 = LBOUND(a, 3), UBOUND(a, 3)
          s = add_line(s, serialize(a(k1, k2, k3)))
        END DO
      END DO
    END DO
    IF (LEN(s) > 0) s = s(: LEN(s) - 1)
  END FUNCTION integer__2_2s3
  FUNCTION integer__4_2s3(a) RESULT(s)
    INTEGER(KIND = 4), INTENT(IN) :: a(:, :, :)
    CHARACTER(LEN = :), ALLOCATABLE :: s
    INTEGER :: k, k1, k2, k3
    s = ""
    s = add_line(s, "# entries")
    DO k1 = LBOUND(a, 1), UBOUND(a, 1)
      DO k2 = LBOUND(a, 2), UBOUND(a, 2)
        DO k3 = LBOUND(a, 3), UBOUND(a, 3)
          s = add_line(s, serialize(a(k1, k2, k3)))
        END DO
      END DO
    END DO
    IF (LEN(s) > 0) s = s(: LEN(s) - 1)
  END FUNCTION integer__4_2s3
  FUNCTION integer__8_2s3(a) RESULT(s)
    INTEGER(KIND = 8), INTENT(IN) :: a(:, :, :)
    CHARACTER(LEN = :), ALLOCATABLE :: s
    INTEGER :: k, k1, k2, k3
    s = ""
    s = add_line(s, "# entries")
    DO k1 = LBOUND(a, 1), UBOUND(a, 1)
      DO k2 = LBOUND(a, 2), UBOUND(a, 2)
        DO k3 = LBOUND(a, 3), UBOUND(a, 3)
          s = add_line(s, serialize(a(k1, k2, k3)))
        END DO
      END DO
    END DO
    IF (LEN(s) > 0) s = s(: LEN(s) - 1)
  END FUNCTION integer__8_2s3
  FUNCTION real__4_2s3(a) RESULT(s)
    REAL(KIND = 4), INTENT(IN) :: a(:, :, :)
    CHARACTER(LEN = :), ALLOCATABLE :: s
    INTEGER :: k, k1, k2, k3
    s = ""
    s = add_line(s, "# entries")
    DO k1 = LBOUND(a, 1), UBOUND(a, 1)
      DO k2 = LBOUND(a, 2), UBOUND(a, 2)
        DO k3 = LBOUND(a, 3), UBOUND(a, 3)
          s = add_line(s, serialize(a(k1, k2, k3)))
        END DO
      END DO
    END DO
    IF (LEN(s) > 0) s = s(: LEN(s) - 1)
  END FUNCTION real__4_2s3
  FUNCTION real__8_2s3(a) RESULT(s)
    REAL(KIND = 8), INTENT(IN) :: a(:, :, :)
    CHARACTER(LEN = :), ALLOCATABLE :: s
    INTEGER :: k, k1, k2, k3
    s = ""
    s = add_line(s, "# entries")
    DO k1 = LBOUND(a, 1), UBOUND(a, 1)
      DO k2 = LBOUND(a, 2), UBOUND(a, 2)
        DO k3 = LBOUND(a, 3), UBOUND(a, 3)
          s = add_line(s, serialize(a(k1, k2, k3)))
        END DO
      END DO
    END DO
    IF (LEN(s) > 0) s = s(: LEN(s) - 1)
  END FUNCTION real__8_2s3
  FUNCTION logical_2s4(a) RESULT(s)
    LOGICAL, INTENT(IN) :: a(:, :, :, :)
    CHARACTER(LEN = :), ALLOCATABLE :: s
    INTEGER :: k, k1, k2, k3, k4
    s = ""
    s = add_line(s, "# entries")
    DO k1 = LBOUND(a, 1), UBOUND(a, 1)
      DO k2 = LBOUND(a, 2), UBOUND(a, 2)
        DO k3 = LBOUND(a, 3), UBOUND(a, 3)
          DO k4 = LBOUND(a, 4), UBOUND(a, 4)
            s = add_line(s, serialize(a(k1, k2, k3, k4)))
          END DO
        END DO
      END DO
    END DO
    IF (LEN(s) > 0) s = s(: LEN(s) - 1)
  END FUNCTION logical_2s4
  FUNCTION integer__1_2s4(a) RESULT(s)
    INTEGER(KIND = 1), INTENT(IN) :: a(:, :, :, :)
    CHARACTER(LEN = :), ALLOCATABLE :: s
    INTEGER :: k, k1, k2, k3, k4
    s = ""
    s = add_line(s, "# entries")
    DO k1 = LBOUND(a, 1), UBOUND(a, 1)
      DO k2 = LBOUND(a, 2), UBOUND(a, 2)
        DO k3 = LBOUND(a, 3), UBOUND(a, 3)
          DO k4 = LBOUND(a, 4), UBOUND(a, 4)
            s = add_line(s, serialize(a(k1, k2, k3, k4)))
          END DO
        END DO
      END DO
    END DO
    IF (LEN(s) > 0) s = s(: LEN(s) - 1)
  END FUNCTION integer__1_2s4
  FUNCTION integer__2_2s4(a) RESULT(s)
    INTEGER(KIND = 2), INTENT(IN) :: a(:, :, :, :)
    CHARACTER(LEN = :), ALLOCATABLE :: s
    INTEGER :: k, k1, k2, k3, k4
    s = ""
    s = add_line(s, "# entries")
    DO k1 = LBOUND(a, 1), UBOUND(a, 1)
      DO k2 = LBOUND(a, 2), UBOUND(a, 2)
        DO k3 = LBOUND(a, 3), UBOUND(a, 3)
          DO k4 = LBOUND(a, 4), UBOUND(a, 4)
            s = add_line(s, serialize(a(k1, k2, k3, k4)))
          END DO
        END DO
      END DO
    END DO
    IF (LEN(s) > 0) s = s(: LEN(s) - 1)
  END FUNCTION integer__2_2s4
  FUNCTION integer__4_2s4(a) RESULT(s)
    INTEGER(KIND = 4), INTENT(IN) :: a(:, :, :, :)
    CHARACTER(LEN = :), ALLOCATABLE :: s
    INTEGER :: k, k1, k2, k3, k4
    s = ""
    s = add_line(s, "# entries")
    DO k1 = LBOUND(a, 1), UBOUND(a, 1)
      DO k2 = LBOUND(a, 2), UBOUND(a, 2)
        DO k3 = LBOUND(a, 3), UBOUND(a, 3)
          DO k4 = LBOUND(a, 4), UBOUND(a, 4)
            s = add_line(s, serialize(a(k1, k2, k3, k4)))
          END DO
        END DO
      END DO
    END DO
    IF (LEN(s) > 0) s = s(: LEN(s) - 1)
  END FUNCTION integer__4_2s4
  FUNCTION integer__8_2s4(a) RESULT(s)
    INTEGER(KIND = 8), INTENT(IN) :: a(:, :, :, :)
    CHARACTER(LEN = :), ALLOCATABLE :: s
    INTEGER :: k, k1, k2, k3, k4
    s = ""
    s = add_line(s, "# entries")
    DO k1 = LBOUND(a, 1), UBOUND(a, 1)
      DO k2 = LBOUND(a, 2), UBOUND(a, 2)
        DO k3 = LBOUND(a, 3), UBOUND(a, 3)
          DO k4 = LBOUND(a, 4), UBOUND(a, 4)
            s = add_line(s, serialize(a(k1, k2, k3, k4)))
          END DO
        END DO
      END DO
    END DO
    IF (LEN(s) > 0) s = s(: LEN(s) - 1)
  END FUNCTION integer__8_2s4
  FUNCTION real__4_2s4(a) RESULT(s)
    REAL(KIND = 4), INTENT(IN) :: a(:, :, :, :)
    CHARACTER(LEN = :), ALLOCATABLE :: s
    INTEGER :: k, k1, k2, k3, k4
    s = ""
    s = add_line(s, "# entries")
    DO k1 = LBOUND(a, 1), UBOUND(a, 1)
      DO k2 = LBOUND(a, 2), UBOUND(a, 2)
        DO k3 = LBOUND(a, 3), UBOUND(a, 3)
          DO k4 = LBOUND(a, 4), UBOUND(a, 4)
            s = add_line(s, serialize(a(k1, k2, k3, k4)))
          END DO
        END DO
      END DO
    END DO
    IF (LEN(s) > 0) s = s(: LEN(s) - 1)
  END FUNCTION real__4_2s4
  FUNCTION real__8_2s4(a) RESULT(s)
    REAL(KIND = 8), INTENT(IN) :: a(:, :, :, :)
    CHARACTER(LEN = :), ALLOCATABLE :: s
    INTEGER :: k, k1, k2, k3, k4
    s = ""
    s = add_line(s, "# entries")
    DO k1 = LBOUND(a, 1), UBOUND(a, 1)
      DO k2 = LBOUND(a, 2), UBOUND(a, 2)
        DO k3 = LBOUND(a, 3), UBOUND(a, 3)
          DO k4 = LBOUND(a, 4), UBOUND(a, 4)
            s = add_line(s, serialize(a(k1, k2, k3, k4)))
          END DO
        END DO
      END DO
    END DO
    IF (LEN(s) > 0) s = s(: LEN(s) - 1)
  END FUNCTION real__8_2s4
  FUNCTION logical_2s(x) RESULT(s)
    LOGICAL, INTENT(IN) :: x
    INTEGER :: y
    CHARACTER(LEN = :), ALLOCATABLE :: s
    ALLOCATE(CHARACTER(LEN = 50)::s)
    y = x
    WRITE(s, '(g0)') y
    s = TRIM(s)
  END FUNCTION logical_2s
  FUNCTION integer1_2s(x) RESULT(s)
    INTEGER(KIND = 1), INTENT(IN) :: x
    CHARACTER(LEN = :), ALLOCATABLE :: s
    ALLOCATE(CHARACTER(LEN = 50)::s)
    WRITE(s, '(g0)') x
    s = TRIM(s)
  END FUNCTION integer1_2s
  FUNCTION integer2_2s(x) RESULT(s)
    INTEGER(KIND = 2), INTENT(IN) :: x
    CHARACTER(LEN = :), ALLOCATABLE :: s
    ALLOCATE(CHARACTER(LEN = 50)::s)
    WRITE(s, '(g0)') x
    s = TRIM(s)
  END FUNCTION integer2_2s
  FUNCTION integer4_2s(x) RESULT(s)
    INTEGER(KIND = 4), INTENT(IN) :: x
    CHARACTER(LEN = :), ALLOCATABLE :: s
    ALLOCATE(CHARACTER(LEN = 50)::s)
    WRITE(s, '(g0)') x
    s = TRIM(s)
  END FUNCTION integer4_2s
  FUNCTION integer8_2s(x) RESULT(s)
    INTEGER(KIND = 8), INTENT(IN) :: x
    CHARACTER(LEN = :), ALLOCATABLE :: s
    ALLOCATE(CHARACTER(LEN = 50)::s)
    WRITE(s, '(g0)') x
    s = TRIM(s)
  END FUNCTION integer8_2s
  FUNCTION real4_2s(x) RESULT(s)
    REAL(KIND = 4), INTENT(IN) :: x
    CHARACTER(LEN = :), ALLOCATABLE :: s
    ALLOCATE(CHARACTER(LEN = 50)::s)
    WRITE(s, '(g0)') x
    s = TRIM(s)
  END FUNCTION real4_2s
  FUNCTION real8_2s(x) RESULT(s)
    REAL(KIND = 8), INTENT(IN) :: x
    CHARACTER(LEN = :), ALLOCATABLE :: s
    ALLOCATE(CHARACTER(LEN = 50)::s)
    WRITE(s, '(g0)') x
    s = TRIM(s)
  END FUNCTION real8_2s
END MODULE serdesl