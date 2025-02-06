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

#include "icon_definitions.inc"
MODULE mo_communication_factory
  USE mo_exception, ONLY: finish
  USE mo_communication_types, ONLY: t_comm_pattern, t_comm_pattern_collection, &
    & t_p_comm_pattern, xfer_list
  USE mo_communication_orig,   ONLY: t_comm_pattern_orig, &
    &                                t_comm_pattern_collection_orig
#ifdef HAVE_YAXT
  USE mo_communication_yaxt,   ONLY: t_comm_pattern_yaxt, &
    &                                t_comm_pattern_collection_yaxt
#endif
  USE mo_decomposition_tools,  ONLY: t_glb2loc_index_lookup
  USE mo_parallel_config, ONLY: comm_pattern_type_orig, &
    comm_pattern_type_yaxt, default_comm_pattern_type
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: setup_comm_pattern, setup_comm_pattern2
  PUBLIC :: setup_comm_pattern_collection
  CHARACTER(len=*), PARAMETER :: modname &
    = 'mo_communication_factory'
CONTAINS
    !>
  !! Sets up a communication pattern for exchanging data.
  !!
  !! Note: This setup routine works only for the trivial communication
  !!       patterns in sequential runs.
  !!
  !! dst_n_points     Total number of points in the RECEIVER array,
  !!                  not every point is necessarily set during exchange
  !!                  (see owner!)
  !!
  !! dst_owner        Owner PE number of every point in the RECEIVER array,
  !!                  if owner(.) == -1, this point will not be set during exchange.
  !!                  If owner(.) == p_pe, this point will be exchanged,
  !!                  this is necessary if sender and receiver arrays are
  !!                  different (e.g. feedback, gather, scatter)
  !!
  !! dst_global_index Global index of of every point in the RECEIVER array
  !!                  There may be more than 1 point with the same global index,
  !!                  in this case the point is exchanged only once and
  !!                  locally distributed.
  !!                  - If this argument is not present, we assume global_index=1,2.3,...
  !! inplace          In-place data exchanges are allowed (source and destination
  !!                  arrays can be identically for the exchange)
  !!                  - if inplace == true, the user guarantees that
  !!                    (src_n_points == dst_n_points) and that points, which will
  !!                    have to be sent to other processes, are disjunct from the
  !!                    points that will have to be received
  !!                  - in case the user only provides the receive array to an
  !!                    exchange call and not send array, the exchange will be
  !!                    faster if inplace == true
  !!
  !! send_decomp_info domain decomposition information for the SENDER array
  !!
  SUBROUTINE setup_comm_pattern(dst_n_points, dst_owner, dst_global_index, &
                                send_glb2loc_index, src_n_points, src_owner, &
                                src_global_index, p_pat, inplace, comm, &
                                opt_comm_type)
    INTEGER, INTENT(IN) :: dst_n_points        ! Total number of points
    INTEGER, INTENT(IN) :: dst_owner(:)        ! Owner of every point
    INTEGER, INTENT(IN) :: dst_global_index(:) ! Global index of every point
    TYPE(t_glb2loc_index_lookup), INTENT(IN) :: send_glb2loc_index
                                               ! global to local index
                                               ! lookup information
                                               ! of the SENDER array
    INTEGER, INTENT(IN) :: src_n_points        ! Total number of points
    INTEGER, INTENT(IN) :: src_owner(:)        ! Owner of every point
    INTEGER, INTENT(IN) :: src_global_index(:) ! Global index of every point


    CLASS(t_comm_pattern), POINTER, INTENT(INOUT) :: p_pat

    LOGICAL, OPTIONAL, INTENT(IN) :: inplace

    INTEGER, OPTIONAL, INTENT(IN) :: comm, opt_comm_type

    INTEGER :: comm_type
    CHARACTER(len=*), PARAMETER :: &
      routine = modname//"::setup_comm_pattern"

    !-----------------------------------------------------------------------

    IF (PRESENT(opt_comm_type)) THEN
       comm_type = opt_comm_type
    ELSE
       comm_type = default_comm_pattern_type
    END IF

    SELECT CASE (comm_type)
      CASE (comm_pattern_type_orig)
        ALLOCATE(t_comm_pattern_orig::p_pat)
      CASE (comm_pattern_type_yaxt)
#ifdef HAVE_YAXT
        ALLOCATE(t_comm_pattern_yaxt::p_pat)
#else
      CALL finish(routine, &
        "comm_pattern_type_yaxt has been selected, but ICON was built without YAXT")
#endif
      CASE DEFAULT
        CALL finish(routine, "Invalid comm_type!")
    END SELECT

    CALL p_pat%setup(dst_n_points, dst_owner, dst_global_index, &
                     send_glb2loc_index, src_n_points, src_owner, &
                     src_global_index, inplace, comm)

  END SUBROUTINE setup_comm_pattern

  
  !-------------------------------------------------------------------------


  SUBROUTINE setup_comm_pattern2(p_pat, comm, recv_msg, send_msg, &
       glb2loc_index_recv, glb2loc_index_send, opt_comm_type)
    CLASS(t_comm_pattern), POINTER, INTENT(INOUT) :: p_pat
    INTEGER, INTENT(in) :: comm
    TYPE(xfer_list), INTENT(in), CONTIGUOUS :: recv_msg(:), send_msg(:)
    TYPE(t_glb2loc_index_lookup), INTENT(IN) :: glb2loc_index_recv, &
         glb2loc_index_send

   INTEGER, OPTIONAL, INTENT(IN) :: opt_comm_type

   CHARACTER(len=*), PARAMETER :: routine = modname//"::setup_comm_pattern2"
   INTEGER :: comm_type

!-----------------------------------------------------------------------

   IF (PRESENT(opt_comm_type)) THEN
      comm_type = opt_comm_type
   ELSE
      comm_type = default_comm_pattern_type
   END IF

   SELECT CASE (comm_type)
     CASE (comm_pattern_type_orig)
      ALLOCATE(t_comm_pattern_orig::p_pat)
     CASE (comm_pattern_type_yaxt)
#ifdef HAVE_YAXT
       ALLOCATE(t_comm_pattern_yaxt::p_pat)
#else
       CALL finish(routine, &
         "comm_pattern_type_yaxt has been selected, but ICON was built without YAXT")
#endif
     CASE DEFAULT
       CALL finish(routine, "Invalid comm_type!")
   END SELECT

   CALL p_pat%setup2(comm, recv_msg, send_msg, &
     &               glb2loc_index_recv, glb2loc_index_send)

  END SUBROUTINE setup_comm_pattern2

  SUBROUTINE setup_comm_pattern_collection(patterns, pattern_collection)

    TYPE(t_p_comm_pattern), INTENT(IN) :: patterns(:)
    CLASS(t_comm_pattern_collection), POINTER :: pattern_collection

    CLASS(t_comm_pattern), POINTER :: first_pattern
    INTEGER :: i
    CHARACTER(len=*), PARAMETER :: &
      routine = modname//"::setup_comm_pattern_collection"
    IF (SIZE(patterns) == 0) THEN
      CALL finish(routine, "SIZE(patterns) == 0")
    END IF

    first_pattern => patterns(1)%p

    DO i = 2, SIZE(patterns)
      IF (.NOT.(SAME_TYPE_AS(first_pattern, patterns(i)%p))) THEN
        CALL finish(routine, "different pattern types")
      END IF
    END DO

    SELECT TYPE(first_pattern)
      TYPE IS(t_comm_pattern_orig)
        ALLOCATE(t_comm_pattern_collection_orig::pattern_collection)
#ifdef HAVE_YAXT
      TYPE IS(t_comm_pattern_yaxt)
        ALLOCATE(t_comm_pattern_collection_yaxt::pattern_collection)
#endif
      CLASS DEFAULT
        CALL finish(routine, "unknown comm pattern class")
    END SELECT

    CALL pattern_collection%setup(patterns)

  END SUBROUTINE setup_comm_pattern_collection


END MODULE mo_communication_factory
!
! Local Variables:
! f90-continuation-indent: 2
! End:
!
