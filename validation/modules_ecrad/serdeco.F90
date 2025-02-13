MODULE serdeco
  IMPLICIT NONE
  INTERFACE serialize
    MODULE PROCEDURE :: character_2s
    MODULE PROCEDURE cloud_type_2s, cloud_optics_type_2s, config_type_2s, thermodynamics_type_2s, logical_2s1, integer__1_2s1, integer__2_2s1, integer__4_2s1, integer__8_2s1, real__4_2s1, real__8_2s1, logical_2s2, integer__1_2s2, integer__2_2s2, integer__4_2s2, integer__8_2s2, real__4_2s2, real__8_2s2, logical_2s3, integer__1_2s3, integer__2_2s3, integer__4_2s3, integer__8_2s3, real__4_2s3, real__8_2s3, logical_2s4, integer__1_2s4, integer__2_2s4, integer__4_2s4, integer__8_2s4, real__4_2s4, real__8_2s4, logical_2s, integer1_2s, integer2_2s, integer4_2s, integer8_2s, real4_2s, real8_2s
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
    s = add_line(s, '# mixing_ratio')
    s = add_line(s, '# alloc')
    s = add_line(s, serialize(ALLOCATED(x % mixing_ratio)))
    IF (ALLOCATED(x % mixing_ratio)) THEN
      s = add_line(s, "# rank")
      s = add_line(s, serialize(3))
      s = add_line(s, "# size")
      DO kmeta = 1, 3
        s = add_line(s, serialize(SIZE(x % mixing_ratio, kmeta)))
      END DO
      s = add_line(s, "# lbound")
      DO kmeta = 1, 3
        s = add_line(s, serialize(LBOUND(x % mixing_ratio, kmeta)))
      END DO
      s = add_line(s, serialize(x % mixing_ratio))
    END IF
    s = add_line(s, '# effective_radius')
    s = add_line(s, '# alloc')
    s = add_line(s, serialize(ALLOCATED(x % effective_radius)))
    IF (ALLOCATED(x % effective_radius)) THEN
      s = add_line(s, "# rank")
      s = add_line(s, serialize(3))
      s = add_line(s, "# size")
      DO kmeta = 1, 3
        s = add_line(s, serialize(SIZE(x % effective_radius, kmeta)))
      END DO
      s = add_line(s, "# lbound")
      DO kmeta = 1, 3
        s = add_line(s, serialize(LBOUND(x % effective_radius, kmeta)))
      END DO
      s = add_line(s, serialize(x % effective_radius))
    END IF
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
    IF (LEN(s) > 0) s = s(: LEN(s) - 1)
  END FUNCTION cloud_type_2s
  FUNCTION cloud_optics_type_2s(x) RESULT(s)
    USE radiation_cloud_optics_data, ONLY: cloud_optics_type
    TYPE(cloud_optics_type), TARGET, INTENT(IN) :: x
    CHARACTER(LEN = :), ALLOCATABLE :: s
    INTEGER :: kmeta, kmeta_0, kmeta_1, kmeta_2, kmeta_3, kmeta_4, kmeta_5, kmeta_6, kmeta_7, kmeta_8, kmeta_9
    s = ""
    s = add_line(s, '# liq_coeff_lw')
    s = add_line(s, '# alloc')
    s = add_line(s, serialize(ALLOCATED(x % liq_coeff_lw)))
    IF (ALLOCATED(x % liq_coeff_lw)) THEN
      s = add_line(s, "# rank")
      s = add_line(s, serialize(2))
      s = add_line(s, "# size")
      DO kmeta = 1, 2
        s = add_line(s, serialize(SIZE(x % liq_coeff_lw, kmeta)))
      END DO
      s = add_line(s, "# lbound")
      DO kmeta = 1, 2
        s = add_line(s, serialize(LBOUND(x % liq_coeff_lw, kmeta)))
      END DO
      s = add_line(s, serialize(x % liq_coeff_lw))
    END IF
    s = add_line(s, '# liq_coeff_sw')
    s = add_line(s, '# alloc')
    s = add_line(s, serialize(ALLOCATED(x % liq_coeff_sw)))
    IF (ALLOCATED(x % liq_coeff_sw)) THEN
      s = add_line(s, "# rank")
      s = add_line(s, serialize(2))
      s = add_line(s, "# size")
      DO kmeta = 1, 2
        s = add_line(s, serialize(SIZE(x % liq_coeff_sw, kmeta)))
      END DO
      s = add_line(s, "# lbound")
      DO kmeta = 1, 2
        s = add_line(s, serialize(LBOUND(x % liq_coeff_sw, kmeta)))
      END DO
      s = add_line(s, serialize(x % liq_coeff_sw))
    END IF
    s = add_line(s, '# ice_coeff_lw')
    s = add_line(s, '# alloc')
    s = add_line(s, serialize(ALLOCATED(x % ice_coeff_lw)))
    IF (ALLOCATED(x % ice_coeff_lw)) THEN
      s = add_line(s, "# rank")
      s = add_line(s, serialize(2))
      s = add_line(s, "# size")
      DO kmeta = 1, 2
        s = add_line(s, serialize(SIZE(x % ice_coeff_lw, kmeta)))
      END DO
      s = add_line(s, "# lbound")
      DO kmeta = 1, 2
        s = add_line(s, serialize(LBOUND(x % ice_coeff_lw, kmeta)))
      END DO
      s = add_line(s, serialize(x % ice_coeff_lw))
    END IF
    s = add_line(s, '# ice_coeff_sw')
    s = add_line(s, '# alloc')
    s = add_line(s, serialize(ALLOCATED(x % ice_coeff_sw)))
    IF (ALLOCATED(x % ice_coeff_sw)) THEN
      s = add_line(s, "# rank")
      s = add_line(s, serialize(2))
      s = add_line(s, "# size")
      DO kmeta = 1, 2
        s = add_line(s, serialize(SIZE(x % ice_coeff_sw, kmeta)))
      END DO
      s = add_line(s, "# lbound")
      DO kmeta = 1, 2
        s = add_line(s, serialize(LBOUND(x % ice_coeff_sw, kmeta)))
      END DO
      s = add_line(s, serialize(x % ice_coeff_sw))
    END IF
    IF (LEN(s) > 0) s = s(: LEN(s) - 1)
  END FUNCTION cloud_optics_type_2s
  FUNCTION config_type_2s(x) RESULT(s)
    USE radiation_config, ONLY: config_type
    TYPE(config_type), TARGET, INTENT(IN) :: x
    CHARACTER(LEN = :), ALLOCATABLE :: s
    INTEGER :: kmeta, kmeta_0, kmeta_1, kmeta_2, kmeta_3, kmeta_4, kmeta_5, kmeta_6, kmeta_7, kmeta_8, kmeta_9
    s = ""
    s = add_line(s, '# cloud_optics')
    s = add_line(s, serialize(x % cloud_optics))
    IF (LEN(s) > 0) s = s(: LEN(s) - 1)
  END FUNCTION config_type_2s
  FUNCTION thermodynamics_type_2s(x) RESULT(s)
    USE radiation_thermodynamics, ONLY: thermodynamics_type
    TYPE(thermodynamics_type), TARGET, INTENT(IN) :: x
    CHARACTER(LEN = :), ALLOCATABLE :: s
    INTEGER :: kmeta, kmeta_0, kmeta_1, kmeta_2, kmeta_3, kmeta_4, kmeta_5, kmeta_6, kmeta_7, kmeta_8, kmeta_9
    s = ""
    s = add_line(s, '# pressure_hl')
    s = add_line(s, '# alloc')
    s = add_line(s, serialize(ALLOCATED(x % pressure_hl)))
    IF (ALLOCATED(x % pressure_hl)) THEN
      s = add_line(s, "# rank")
      s = add_line(s, serialize(2))
      s = add_line(s, "# size")
      DO kmeta = 1, 2
        s = add_line(s, serialize(SIZE(x % pressure_hl, kmeta)))
      END DO
      s = add_line(s, "# lbound")
      DO kmeta = 1, 2
        s = add_line(s, serialize(LBOUND(x % pressure_hl, kmeta)))
      END DO
      s = add_line(s, serialize(x % pressure_hl))
    END IF
    IF (LEN(s) > 0) s = s(: LEN(s) - 1)
  END FUNCTION thermodynamics_type_2s
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
END MODULE serdeco