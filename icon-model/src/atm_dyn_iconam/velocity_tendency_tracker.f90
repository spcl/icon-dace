MODULE velocity_tendency_tracker
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: velocity_tendencies_c1_count, &
            velocity_tendencies_c2_count, &
            velocity_tendencies_c3_count, &
            velocity_tendencies_c4_count

  INTEGER :: velocity_tendencies_c1_count = 0
  INTEGER :: velocity_tendencies_c2_count = 0
  INTEGER :: velocity_tendencies_c3_count = 0
  INTEGER :: velocity_tendencies_c4_count = 0
END MODULE velocity_tendency_tracker